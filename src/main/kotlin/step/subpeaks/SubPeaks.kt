package step.subpeaks

import model.*
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.*

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

fun <T : GaussianParameters> fit(
    values: List<Double>,
    offset: Int,
    initParameters: (region: Region) -> T,
    optimize: Optimizer<T>,
    parametersToRegion: (parameters: T, offset: Int) -> Region
): List<Fit<T>> {
    // Subtract out background
    val background = values.min()!!
    val valuesWithoutBackground = values.map { it - background }

    val splitValues = splitForFit(valuesWithoutBackground)

    val fits = mutableListOf<Fit<T>>()
    var fitOffset = offset
    for (valuesToFit in splitValues) {
        // Put all on same scale
        val avg = valuesToFit.average()
        val scaledValues = valuesToFit.map { it / avg }

        // Get initial candidates list
        val candidateGaussians = findCandidates(scaledValues, initParameters)

        log.debug { "Candidate gaussians size: ${candidateGaussians.size}" }
        log.debug { "Candidate gaussians: $candidateGaussians" }
        log.debug { "Candidate gaussian regions (sorted): ${candidateGaussians.map { it.region }.sortedBy { it.start }}" }

        /*
         * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
         * optimization is performed with the Levenberg-Marquardt algorithm
         */
        val optimized = optimize(scaledValues.toDoubleArray(), candidateGaussians,
            candidateGaussians.map { it.parameters })

        log.debug { "Optimized Result size: ${optimized.parameters.size}" }
        log.debug { "Optimized Result error: ${optimized.error}" }
        log.debug { "Optimized Result iterations: ${optimized.iterations}" }
        log.debug { "Optimized Result: $optimized" }
        val sqrtAvg = sqrt(avg)
        // Bring parameters back to original scale and add background back in. Also add offset to mean
        optimized.parameters.forEach {
            it.amplitude *= sqrtAvg
            it.mean += fitOffset
        }

        val subPeaks = optimized.parameters
            .sortedBy { it.mean }
            .map { SubPeak(parametersToRegion(it, fitOffset), parametersToScore(it), it) }

        fits += Fit(Region(fitOffset, fitOffset + valuesToFit.size), subPeaks, background, optimized.error)
        fitOffset += valuesToFit.size
    }
    return fits
}

private fun parametersToScore(parameters: GaussianParameters) =
    parameters.amplitude * parameters.amplitude / parameters.stdDev

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

fun splitForFit(values: List<Double>): List<List<Double>> {
    val splitIndex = splitIndexForFit(values) ?: return listOf(values)
    val splitFirst = splitForFit(values.subList(0, splitIndex))
    val splitSecond = splitForFit(values.subList(splitIndex, values.size))
    return splitFirst + splitSecond
}