package step.subpeaks

import model.*
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.*
import step.PDF
import step.Peak
import util.runParallel
import java.util.*

private val log = KotlinLogging.logger {}

interface GaussianParameters {
    var amplitude: Double
    var mean: Double
    var stdDev: Double
}

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

// Minimum number of values in a peak that are worth evaluating
const val MIN_PEAK_VALUES = 10

abstract class Fitter<T : GaussianParameters> (private val name: String, val optimizer: Optimizer<T>) {

    /**
     * Do all sub-peak fits for all peaks on all chromosomes
     */
    fun fitAll(peaks: Map<String, List<Peak>>, pdfs: Map<String, PDF>): Map<String, List<SubPeak<T>>> {
        log.info { "Running $name for ${peaks.size} chromosomes..." }
        val subPeaks = mutableMapOf<String, List<SubPeak<T>>>()
        for ((chr, chrPeaks) in peaks) {
            subPeaks[chr] = fitChrom(chr, chrPeaks, pdfs.getValue(chr))
        }
        log.info { "$name complete!" }
        return subPeaks
    }

    /**
     * Do sub-peak fits for all peaks on a single chromosome
     */
    fun fitChrom(chr: String, peaks: List<Peak>, pdf: PDF): List<SubPeak<T>> {
        (peaks as MutableList).sortByDescending { it.region.end - it.region.start }
        val subPeaks = Collections.synchronizedList(mutableListOf<SubPeak<T>>())
        runParallel("$name on $chr", "peaks", peaks) { peak ->
            val region = peak.region
            val peakValues = (region.start..region.end).map { pdf[it] }
            subPeaks += fitPeak(peakValues, region.start).flatMap { fit ->
                fit.subPeaks
            }
        }
        subPeaks.sortBy { it.region.start }
        return subPeaks
    }

    /**
     * Do sub-peak fits for a single peak
     */
    fun fitPeak(values: List<Double>, offset: Int): List<Fit<T>> {
        if (values.size < MIN_PEAK_VALUES) return listOf()

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
            val candidateRegions = findCandidates(scaledValues, splitValues.size > 1)
            if (candidateRegions.isEmpty()) continue
            val candidateGaussians = candidateGaussians(scaledValues, candidateRegions, ::initParameters)

            /*
             * Perform the actual fit of the curve to a sum of gaussians.size Gaussians
             * optimization is performed with the Levenberg-Marquardt algorithm
             */
            val optimized = optimizer.optimize(scaledValues.toDoubleArray(), candidateGaussians,
                candidateGaussians.map { it.parameters })

            val sqrtAvg = sqrt(avg)
            // Bring parameters back to original scale and add background back in. Also add offset to mean
            optimized.parameters.forEach {
                it.amplitude *= sqrtAvg
                it.mean += fitOffset
            }

            val subPeaks = optimized.parameters
                .sortedBy { it.mean }
                .map { SubPeak(parametersToRegion(it), parametersToScore(it), it) }

            fits += Fit(Region(fitOffset, fitOffset + valuesToFit.size), subPeaks, background, optimized.error)
            fitOffset += valuesToFit.size
        }
        return fits
    }

    abstract fun initParameters(region: Region): T
    abstract fun parametersToRegion(parameters: T): Region
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