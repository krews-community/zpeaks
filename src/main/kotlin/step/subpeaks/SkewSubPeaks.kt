package step.subpeaks

import model.Region
import mu.KotlinLogging
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.special.Erf.erf
import org.apache.commons.math3.util.FastMath.*
import org.apache.commons.math3.util.Precision
import step.PDF
import step.Peak
import util.runParallel
import java.util.*

private val log = KotlinLogging.logger {}

/**
 * Represents a Gaussian distribution with shape (adds skewness).
 * Comparison functions only compare amplitude (used for sorting
 * from most to least important when curve fitting).
 *
 * shape: shape parameter (alpha) capturing skewness.
 * amplitude: suqare root of the amplitude of the distribution.
 * mean: center point of the distribution.
 * stdev: standard deviation of the distribution.
 */
data class SkewGaussianParameters(
    override var amplitude: Double,
    override var mean: Double,
    override var stdDev: Double,
    var shape: Double
) : GaussianParameters

typealias SkewSubPeak = SubPeak<SkewGaussianParameters>

fun runSkewSubPeaks(peaks: Map<String, List<Peak>>, pdfs: Map<String, PDF>): Map<String, List<SkewSubPeak>> {
    log.info { "Running skew sub-Peaks for ${peaks.size} chromosomes..." }
    val subPeaks = mutableMapOf<String, List<SkewSubPeak>>()
    for ((chr, chrPeaks) in peaks) {
        subPeaks[chr] = runChromSkewSubPeaks(chr, chrPeaks, pdfs.getValue(chr))
    }
    log.info { "Skew sub-peaks complete!" }
    return subPeaks
}

fun runChromSkewSubPeaks(chr: String, peaks: List<Peak>, pdf: PDF): List<SkewSubPeak> {
    (peaks as MutableList).sortByDescending { it.region.end - it.region.start }
    val subPeaks = Collections.synchronizedList(mutableListOf<SkewSubPeak>())
    runParallel("Skew Sub-Peaks on $chr", peaks) { peak ->
        val region = peak.region
        val peakValues = (region.start..region.end).map { pdf[it] }
        subPeaks += fitSkew(peakValues, region.start).flatMap { fit ->
            fit.subPeaks
        }
    }
    subPeaks.sortBy { it.region.start }
    return subPeaks
}

fun fitSkew(values: List<Double>, pileUpStart: Int) =
    fit(values, pileUpStart, ::initSkewParameters, ::optimizeSkew, ::skewParametersToRegion)

fun initSkewParameters(region: Region) = SkewGaussianParameters(
    amplitude = 0.0001,
    mean = (region.start + region.end) / 2.0,
    stdDev = (region.end - region.start) / 2.0,
    shape = 0.0
)

val SQRT2 = sqrt(2.0)

/**
 * Calculate curve value at given x coordinate
 */
fun curveValue(x: Int, amplitude: Double, mean: Double, stdDev: Double, shape: Double): Double {
    return amplitude * amplitude / stdDev *
            exp( -(x - mean) * (x - mean) / stdDev / stdDev / 2) *
            (1.0 + erf(shape * (x - mean) / (stdDev * SQRT2)))
}

/**
 * Calculate curve for given parameters (in raw DoubleArray) form
 */
fun calculateSkewCurve(parameters: DoubleArray, curveLength: Int, start: Int = 0): DoubleArray {
    val curve = DoubleArray(curveLength) { 0.0 }
    for (j in 0 until curveLength) {
        for (k in 0 until parameters.size step 4) {
            curve[j] += curveValue(j+start, parameters[k], parameters[k+1], parameters[k+2], parameters[k+3])
        }
    }
    return curve
}

val SQRT_PI_OVER_2 = sqrt(PI / 2.0)

fun calculateSkewJacobian(parameters: DoubleArray, curveLength: Int): Array<DoubleArray> {
    val jacobian = Array(curveLength) { DoubleArray(parameters.size) }
    for (j in 0 until curveLength) {
        for (k in 0 until parameters.size step 4) {
            val a = parameters[k] // amplitude
            val m = parameters[k+1] // mean
            val u = parameters[k+2] // standard deviation
            val s = parameters[k+3] // shape
            val dm = j - m
            val exp = exp(-dm * dm / (u * u * 2.0))
            val erf = 1.0 + erf(s * dm / (u * SQRT2))
            val expx = exp(-s * s * dm * dm / (2 * u * u) - dm * dm / (2 * u * u))
            jacobian[j][k] = 2 * a * exp * erf / u // partial with respect to a
            jacobian[j][k+1] = (a * a * dm * exp * erf / (u * u * u)) -
                    expx * a * a * s / (u * u * SQRT_PI_OVER_2) // partial with respect to m
            jacobian[j][k+2] = (a * a * dm * dm * exp * erf) / (u * u * u * u) -
                    expx * a * a * s * dm / (SQRT_PI_OVER_2 * u * u * u) -
                    a * a * exp * erf / (u * u) // partial with respect to u
            jacobian[j][k+3] = expx * a * a * dm / (SQRT_PI_OVER_2 * u * u) // partial with respect to s
        }
    }
    return jacobian
}

typealias CandidateSkewGaussian = CandidateGaussian<SkewGaussianParameters>

fun optimizeSkew(values: DoubleArray,
                 candidateGaussians: List<CandidateSkewGaussian>,
                 initialGaussians: List<SkewGaussianParameters>):
        OptimizeResults<SkewGaussianParameters> {
    val avg = values.average()

    val optimizer = LevenbergMarquardtOptimizer()
        .withInitialStepBoundFactor(0.05)
        .withCostRelativeTolerance(avg * 1e-4)
        .withParameterRelativeTolerance(avg * 1e-8)
        .withOrthoTolerance(avg * 1e-3)
        .withRankingThreshold(Precision.SAFE_MIN)

    val initialParameters = DoubleArray(initialGaussians.size * 4)
    for ((i, gaussian) in initialGaussians.withIndex()) {
        initialParameters[i*4] = gaussian.amplitude
        initialParameters[i*4+1] = gaussian.mean
        initialParameters[i*4+2] = gaussian.stdDev
        initialParameters[i*4+3] = gaussian.shape
    }

    val problem = LeastSquaresBuilder()
        .maxEvaluations(1000)
        .maxIterations(1000)
        .model(
            { params -> calculateSkewCurve(params, values.size) },
            { params -> calculateSkewJacobian(params, values.size) }
        )
        .parameterValidator { params -> validateSkewParameters(params, candidateGaussians) }
        .start(initialParameters)
        .target(values)
        .build()

    val optimum = optimizer.optimize(problem)
    val rawParameters = optimum.point
    val optimizedParameters = mutableListOf<SkewGaussianParameters>()
    for (j in 0 until optimum.point.dimension step 4) {
        optimizedParameters += SkewGaussianParameters(rawParameters.getEntry(j), rawParameters.getEntry(j+1),
            rawParameters.getEntry(j+2), rawParameters.getEntry(j+3))
    }
    return OptimizeResults(optimizedParameters, optimum.rms, optimum.iterations)
}

const val MAX_SKEW = 10
const val MAX_STD_DEV = 400

fun validateSkewParameters(params: RealVector, candidateGaussians: List<CandidateSkewGaussian>): RealVector {
    val validated = params.copy()
    for (j in 0 until candidateGaussians.size) {
        val amplitude = params.getEntry(j*4)
        val mean = params.getEntry(j*4+1)
        val stdDev = params.getEntry(j*4+2)
        val shape = params.getEntry(j*4+3)

        if (amplitude < 0) {
            validated.setEntry(j*4, candidateGaussians[j].parameters.amplitude)
        }

        val skewMode = skewMode(mean, stdDev, shape)
        if (skewMode < candidateGaussians[j].region.start || skewMode > candidateGaussians[j].region.end) {
            validated.setEntry(j*4+1, candidateGaussians[j].parameters.mean)
            validated.setEntry(j*4+2, candidateGaussians[j].parameters.stdDev)
            validated.setEntry(j*4+3, candidateGaussians[j].parameters.shape)
        }
        if (stdDev <= 0 || stdDev > MAX_STD_DEV) {
            validated.setEntry(j*4+2, candidateGaussians[j].parameters.stdDev)
        }
        if (abs(shape ) > MAX_SKEW) validated.setEntry(j*4+3, 0.0)
    }
    return validated
}

fun skewParametersToRegion(parameters: SkewGaussianParameters): Region {
    val mode = skewMode(parameters.mean, parameters.stdDev, parameters.shape)
    val start = mode - parameters.stdDev
    val stop = mode + parameters.stdDev
    return Region(start.toInt(), stop.toInt())
}

val SQRT_2_OVER_PI = sqrt(2.0 / PI)
val SQRT_2_OVER_PI_CUBED = pow(SQRT_2_OVER_PI, 3)
const val FOUR_MINUS_PI_OVER_2 = (4.0 - PI) / 2.0
const val NEG_2_PI = -2.0 * PI

/**
 * Calculate the mode for a skewed gaussian
 */
fun skewMode(mean: Double, stdDev: Double, shape: Double): Double {
    val delta = shape / sqrt(1.0 + shape * shape)
    val uz = SQRT_2_OVER_PI * delta
    val oz = sqrt(1.0 - uz * uz)
    val skewness = FOUR_MINUS_PI_OVER_2 *
            ((SQRT_2_OVER_PI_CUBED * delta * delta * delta) / pow(1.0 - 2.0 * delta * delta / PI, 1.5))
    return (uz - skewness * oz / 2.0 - sgn(shape) / 2.0 * exp(NEG_2_PI / abs(shape))) *
            stdDev + mean
}

fun sgn(value: Double) = if (value > 0.0) 1 else if (value < 0.0) -1 else 0