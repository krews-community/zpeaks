package step.subpeaks

import model.Region
import mu.KotlinLogging
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.special.Erf.erf
import org.apache.commons.math3.util.FastMath.*

val SQRT2 = sqrt(2.0)
val SQRT_PI_OVER_2 = sqrt(PI / 2.0)
val SQRT_2_OVER_PI = sqrt(2.0 / PI)
val SQRT_2_OVER_PI_CUBED = pow(SQRT_2_OVER_PI, 3)
const val FOUR_MINUS_PI_OVER_2 = (4.0 - PI) / 2.0
const val NEG_2_PI = -2.0 * PI

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

object SkewFitter : Fitter<SkewGaussianParameters>("Fit Skew Sub-Peaks", SkewOptimizer) {
    override fun initParameters(region: Region) = SkewGaussianParameters(
        amplitude = 0.0,
        mean = (region.start + region.end) / 2.0,
        stdDev = (region.end - region.start) / 2.0,
        shape = 0.0
    )

    override fun parametersToRegion(parameters: SkewGaussianParameters): Region {
        val mode = skewMode(parameters.mean, parameters.stdDev, parameters.shape)
        val start = mode - parameters.stdDev
        val stop = mode + parameters.stdDev
        return Region(start.toInt(), stop.toInt())
    }
}

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

object SkewOptimizer : Optimizer<SkewGaussianParameters>() {

    override fun parametersToArray(parameters: List<SkewGaussianParameters>): DoubleArray {
        val paramArray = DoubleArray(parameters.size * 4)
        for ((i, gaussian) in parameters.withIndex()) {
            paramArray[i * 4] = gaussian.amplitude
            paramArray[i * 4 + 1] = gaussian.mean
            paramArray[i * 4 + 2] = gaussian.stdDev
            paramArray[i * 4 + 3] = gaussian.shape
        }
        return paramArray
    }

    override fun arrayToParameters(array: DoubleArray): List<SkewGaussianParameters> {
        val parameters = mutableListOf<SkewGaussianParameters>()
        for (j in 0 until array.size step 4) {
            parameters += SkewGaussianParameters(array[j], array[j + 1], array[j + 2], array[j + 3])
        }
        return parameters
    }

    override fun calculateCurveValue(x: Int, gaussian: SkewGaussianParameters) = with(gaussian) {
        amplitude * amplitude / stdDev *
                exp(-(x - mean) * (x - mean) / stdDev / stdDev / 2) *
                (1.0 + erf(shape * (x - mean) / (stdDev * SQRT2)))
    }

    override fun calculateJacobian(parameters: DoubleArray, curveLength: Int): Array<DoubleArray> {
        val jacobian = Array(curveLength) { DoubleArray(parameters.size) }
        for (j in 0 until curveLength) {
            for (k in 0 until parameters.size step 4) {
                val a = parameters[k] // amplitude
                val m = parameters[k + 1] // mean
                val u = parameters[k + 2] // standard deviation
                val s = parameters[k + 3] // shape
                val dm = j - m
                val exp = exp(-dm * dm / (u * u * 2.0))
                val erf = 1.0 + erf(s * dm / (u * SQRT2))
                val expx = exp(-s * s * dm * dm / (2 * u * u) - dm * dm / (2 * u * u))
                jacobian[j][k] = 2 * a * exp * erf / u // partial with respect to a
                jacobian[j][k + 1] = (a * a * dm * exp * erf / (u * u * u)) -
                        expx * a * a * s / (u * u * SQRT_PI_OVER_2) // partial with respect to m
                jacobian[j][k + 2] = (a * a * dm * dm * exp * erf) / (u * u * u * u) -
                        expx * a * a * s * dm / (SQRT_PI_OVER_2 * u * u * u) -
                        a * a * exp * erf / (u * u) // partial with respect to u
                jacobian[j][k + 3] = expx * a * a * dm / (SQRT_PI_OVER_2 * u * u) // partial with respect to s
            }
        }
        return jacobian
    }

    private const val MAX_SKEW = 10
    private const val MAX_STD_DEV = 400

    override fun validateParameters(
        parameters: RealVector,
        candidateGaussians: List<CandidateGaussian<SkewGaussianParameters>>
    ): RealVector {
        val validated = parameters.copy()
        for (j in 0 until candidateGaussians.size) {
            val amplitude = parameters.getEntry(j * 4)
            val mean = parameters.getEntry(j * 4 + 1)
            val stdDev = parameters.getEntry(j * 4 + 2)
            val shape = parameters.getEntry(j * 4 + 3)

            if (amplitude <= 0) {
                validated.setEntry(j * 4, candidateGaussians[j].parameters.amplitude)
            }

            val skewMode = skewMode(mean, stdDev, shape)
            if (skewMode < candidateGaussians[j].region.start || skewMode > candidateGaussians[j].region.end) {
                validated.setEntry(j * 4 + 1, candidateGaussians[j].parameters.mean)
                validated.setEntry(j * 4 + 2, candidateGaussians[j].parameters.stdDev)
                validated.setEntry(j * 4 + 3, candidateGaussians[j].parameters.shape)
            }
            if (stdDev <= 0 || stdDev > MAX_STD_DEV) {
                validated.setEntry(j * 4 + 2, candidateGaussians[j].parameters.stdDev)
            }
            if (abs(shape) > MAX_SKEW) validated.setEntry(j * 4 + 3, 0.0)
        }
        return validated
    }

}
