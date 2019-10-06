package step.subpeaks

import model.*
import mu.KotlinLogging
import step.*
import util.*
import java.util.*
import kotlin.math.pow

private val log = KotlinLogging.logger {}

data class ReplicatedFit<T : GaussianParameters>(
    val region: Region,
    val subPeaks: List<ReplicatedSubPeak<T>>,
    val background: Double,
    val error: Double
)

data class ReplicatedSubPeak<T : GaussianParameters>(
    val region: Region,
    val score: Double,
    val gaussianParameters: T,
    val replicationScore: Double
)

abstract class ReplicatedFitter<T : GaussianParameters> (private val name: String, override val optimizer: Optimizer<T>)
    : Fitter<T>(name, optimizer) {

    /**
     * Do sub-peak fits for all peaks on a single chromosome
     */
    fun fit(chr: String, peaks: List<Region>, pdf: List<PDF>): List<ReplicatedSubPeak<T>> {
        (peaks as MutableList).sortByDescending { it.end - it.start }
        val subPeaks = Collections.synchronizedList(mutableListOf<ReplicatedSubPeak<T>>())
        runParallel("$name on $chr", "peaks", peaks) { peak ->
            val peakValues = pdf.map { item -> (peak.start..peak.end).map { item[it].toDouble() } }
            try {
                subPeaks += fitReplicatedPeaks(peakValues, peak.start).flatMap { fit ->
                    fit.subPeaks
                }
            } catch(e: Exception) {
                println("When fitting peak $peak, caught exception: $e")
            }
        }
        subPeaks.sortBy { it.region.start }
        return subPeaks
    }

    private fun contribution(values: List<List<Double>>, fitResult: Fit<T>): ReplicatedFit<T> {
        val parameters = fitResult.subPeaks.sortedByDescending { it.score }
        val currentValues = optimizer.calculateCurve(parameters.map { it.gaussianParameters }, values.size)
        return ReplicatedFit<T>(
            region = fitResult.region, background = fitResult.background, error = fitResult.error,
            subPeaks = fitResult.subPeaks.map {
                val thisCurve = optimizer.calculateCurve(listOf(it.gaussianParameters), values.size)
                val thisRemoved = currentValues.mapIndexed { index, it ->
                    it - thisCurve[index]
                }
                ReplicatedSubPeak<T>(
                    region = it.region, score = it.score, gaussianParameters = it.gaussianParameters,
                    replicationScore = values.map { replicate ->
                        thisRemoved.foldIndexed(0.0, { index, acc, it ->
                            acc + (it - replicate[index]).pow(2.0)
                        }) - currentValues.foldIndexed(0.0, { index, acc, it ->
                            acc + (it - replicate[index]).pow(2.0)
                        })
                    }.average()
                )
            }
        )
    }

    /**
     * Do sub-peak fits for a single peak
     */
    private fun fitReplicatedPeaks(values: List<List<Double>>, offset: Int): List<ReplicatedFit<T>> {

        // Make sure the list from each replicate has the same length, and that it is long enough
        val length = values[0].size
        values.forEach { assert(it.size == length) }
        if (length < MIN_PEAK_VALUES) return listOf()

        // Take the average and fit
        val averageValues: List<Double> = values.fold(values[0], { acc, it ->
            it.mapIndexed { index, element -> acc[index] + element }
        })
        return fitPeak(averageValues, offset).map { contribution(values, it) }

    }

}

private fun parametersToScore(parameters: GaussianParameters) =
    parameters.amplitude * parameters.amplitude / parameters.stdDev

object StandardReplicatedFitter : ReplicatedFitter<StandardGaussianParameters>("Fit Standard Sub-Peaks", StandardOptimizer) {
    override fun initParameters(region: Region) = StandardGaussianParameters(
        amplitude = 0.0,
        mean = (region.start + region.end) / 2.0,
        stdDev = (region.end - region.start) / 2.0
    )

    override fun parametersToRegion(parameters: StandardGaussianParameters) = Region(
        start = (parameters.mean - parameters.stdDev).toInt(),
        end = (parameters.mean + parameters.stdDev).toInt()
    )
}

object SkewReplicatedFitter : ReplicatedFitter<SkewGaussianParameters>("Fit Skew Sub-Peaks", SkewOptimizer) {
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
