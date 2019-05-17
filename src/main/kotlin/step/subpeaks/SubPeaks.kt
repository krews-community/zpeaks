package step.subpeaks

import model.*
import mu.KotlinLogging
import org.apache.commons.math3.linear.*
import org.apache.commons.math3.linear.MatrixUtils.*
import org.apache.commons.math3.util.FastMath.*
import util.*

private val log = KotlinLogging.logger {}

interface GaussianParameters {
    var amplitude: Double
    var mean: Double
    var stdDev: Double
}

data class CandidateGaussian<T : GaussianParameters>(
    val region: Region,
    val parameters: T
)

data class Fit<T : GaussianParameters>(
    val region: Region,
    val subPeaks: List<SubPeak<T>>,
    val background: Double,
    val error: Double
)

data class SubPeak<T : GaussianParameters>(
    val region: Region,
    val score: Double,
    val gaussianParameters: T
)

data class OptimizeResults<T : GaussianParameters>(
    val parameters: List<T>,
    val error: Double,
    val iterations: Int
)

typealias Optimizer<T> = (values: DoubleArray, candidateGaussians: List<CandidateGaussian<T>>, gaussians: List<T>) -> OptimizeResults<T>

const val SUB_PEAKS_HARD_MAX = 5000
const val SUB_PEAKS_SOFT_MAX = 2000
const val SUB_PEAKS_SOFT_MAX_RATIO = 0.05

fun splitIndexForFit(values: List<Double>): Int? {
    if (values.size <= SUB_PEAKS_SOFT_MAX) return null

    // look only within a buffer to make sure we shrink the sizes by a reasonable amount
    val bufferSize = values.size / 5
    val minIndex = values.subList(bufferSize, values.size - bufferSize)
        .withIndex().minBy { it.value }!!.index + bufferSize

    if (values.size <= SUB_PEAKS_HARD_MAX) {
        val minMaxRatio = values[minIndex] / values.max()!!
        if (minMaxRatio > SUB_PEAKS_SOFT_MAX_RATIO) return null
    }

    // Return two halves, split at min index.
    return minIndex
}

fun <T : GaussianParameters> fit(
    values: List<Double>,
    offset: Int,
    initParameters: (region: Region) -> T,
    optimize: Optimizer<T>,
    parametersToRegion: (parameters: T, offset: Int) -> Region
): List<Fit<T>> {
    // Subtract out background
    val background = max(values.first(), values.last())
    val valuesWithoutBackground = values.map { it - background }

    val splitIndex = splitIndexForFit(valuesWithoutBackground)
    if (splitIndex != null) {
        val splitFirst = values.subList(0, splitIndex)
        val splitSecond = values.subList(splitIndex, values.size)
        log.info { "Splitting: ${values.size} to ${splitFirst.size} and ${splitSecond.size} for offset $offset" }
        val firstFits = fit(splitFirst, offset, initParameters, optimize, parametersToRegion)
        val secondFits = fit(splitSecond, offset + splitFirst.size,
            initParameters, optimize, parametersToRegion)
        return firstFits + secondFits
    }

    // In order for the findCandidates algorithm to work properly, the start and end values
    // should be similar, but when we split we may make them lob-sided. fitRegion cuts one side to correct this.
    val fitRegionStart =
        if (values.first() == background) 0
        else values.withIndex().first { it.value >= background }.index
    val fitRegionEnd =
        if (values.last() == background) values.size
        else values.withIndex().last { it.value >= background }.index
    val fitRegion = Region(fitRegionStart, fitRegionEnd)
    // This is the fit region offset by the start. It represents a region with absolute positions in the chromosome.
    val offsetFitRegion = Region(fitRegion.start + offset, fitRegion.end + offset)
    val fitValues = valuesWithoutBackground.subList(fitRegion.start, fitRegion.end)

    // Put all on same scale
    val avg = fitValues.average()
    val scaledValues = fitValues.map { it / avg }

    // Get initial candidates list
    val candidates = findCandidates(scaledValues)

    /*
     * here we compute guesses for the Gaussians' parameters based on zero crossing locations
     * the means and standard deviations are computed with simple equations
     * the amplitudes are fit to the data using least squares after the means and standard deviations are computed
     */
    val candidateGaussians = candidates
        .map { CandidateGaussian(it, initParameters(it)) } as MutableList

    /*
     * Solve the system of linear equations to get the amplitude guesses
     */
    val lim = scaledValues.size
    val m = Array(candidateGaussians.size) { DoubleArray(candidateGaussians.size) }
    val d = DoubleArray(candidateGaussians.size)
    val distributions = candidateGaussians
        .map { gaussianDistribution(1.0, it.parameters.mean, it.parameters.stdDev, lim) }

    // Set coefficients: sum all bp x of (sum all gaussians j of (kth gaussian(x) * jth gaussian(x)))
    for (j in 0 until candidateGaussians.size) {
        for (k in 0 until candidateGaussians.size) {
            var value = 0.0
            for (n in 0 until lim) {
                value += distributions[k][n] * distributions[j][n]
            }
            m[k][j] = value
        }
    }

    // Set values: sum all bp x of (curve value x * jth gaussian(x))
    for (j in 0 until candidateGaussians.size) {
        var value = 0.0
        for (n in 0 until lim) {
            value += distributions[j][n] * scaledValues[n]
        }
        d[j] = value
    }

    val r = LUDecomposition(createRealMatrix(m)).solver.solve(createRealVector(d))

    for (j in 0 until candidateGaussians.size) {
        val s = if (r.getEntry(j) < 0) -1.0 else 1.0
        // Set Amplitude
        candidateGaussians[j].parameters.amplitude =
            sqrt(abs(r.getEntry(j)) * candidateGaussians[j].parameters.stdDev) * s
    }
    // Sort gaussians by amplitude
    candidateGaussians.sortByDescending { it.parameters.amplitude * it.parameters.amplitude / it.parameters.stdDev }

    log.debug { "Candidate gaussians size: ${candidateGaussians.size}" }
    log.debug { "Candidate gaussians: $candidateGaussians" }
    log.debug { "Candidate gaussian regions (sorted): ${candidateGaussians.map { it.region }.sortedBy { it.start }}" }

    /*
     * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
     * optimization is performed with the Levenberg-Marquardt algorithm
     */
    var bestOptimized: OptimizeResults<T>? = null
    var previousOptimized: OptimizeResults<T>? = null
    for (j in candidateGaussians.size / 2 until candidateGaussians.size) {
        val currentGuess =
            candidateGaussians.subList(0, j + 1).map { it.parameters }
        //    (previousOptimized?.parameters ?: candidateGaussians.subList(0, j + 1).map { it.parameters }) + candidateGaussians[j].parameters

        val optimizeResults = optimize(scaledValues.toDoubleArray(), candidateGaussians.subList(0, j + 1), currentGuess)
        log.debug { "Attempt #$j. Current guess size = ${currentGuess.size}. Error = ${optimizeResults.error}" }

        /*previousOptimized = optimizeResults
        if (optimizeResults.error < 0.05) {
            log.debug { "Answer has acceptable error < 0.05!" }
            bestOptimized = optimizeResults
            break
        }*/
        if (bestOptimized == null || optimizeResults.error < bestOptimized.error) {
            bestOptimized = optimizeResults
        }
    }

    log.debug { "Optimized Result size: ${bestOptimized!!.parameters.size}" }
    log.debug { "Optimized Result error: ${bestOptimized!!.error}" }
    log.debug { "Optimized Result iterations: ${bestOptimized!!.iterations}" }
    log.debug { "Optimized Result: $bestOptimized" }
    val sqrtAvg = sqrt(avg)
    // Bring parameters back to original scale and add background back in. Also add offset to mean
    bestOptimized!!.parameters.forEach {
        it.amplitude *= sqrtAvg
        it.mean += offsetFitRegion.start
    }

    val subPeaks = bestOptimized.parameters
        .sortedBy { it.mean }
        .map { SubPeak(parametersToRegion(it, offsetFitRegion.start), parametersToScore(it), it) }

    return listOf(Fit(offsetFitRegion, subPeaks, background, bestOptimized.error))
}

fun parametersToScore(parameters: GaussianParameters): Double = parameters.amplitude.pow(2) / parameters.stdDev