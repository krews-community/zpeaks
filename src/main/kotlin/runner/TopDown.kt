package runner

import model.*
import step.*
import mu.KotlinLogging

private val log = KotlinLogging.logger {}

class TopDownZRunner(runConfig: ZRunConfig) : ZRunner("Top-Down", runConfig) {

    private var allCurrentChrPeaks: List<List<Region>>? = null

    override fun pileUp(chr: String, length: Int, onRange: IntRange?, subsetSize: Int?): PileUp = with(runConfig) {
        log.info { "Creating peaks-only aggregate track..." }
        // Array containing all aggregated data
        val aggregateChrData = DoubleArray(length)
        // Array containing the number of sources the data came from at each BP
        val aggregateSources = IntArray(length)
        val allPeaks = mutableListOf<List<Region>>()
        allCurrentChrPeaks = allPeaks

        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, length, pileUpInput.options)
            val pdf = pdf(pileUp, smoothing, normalizePDF, onRange, subsetSize)
            val peaks = callPeaks(pdf, threshold)
            allPeaks += peaks
            // Add raw pile-up data to aggregate from peak regions only
            for (peak in peaks) {
                for (bp in peak.start .. peak.end) {
                    aggregateChrData[bp] += pileUp.scaledValue(bp)
                    aggregateSources[bp]++
                }
            }
        }

        // Scale the aggregate peaks data based on the number of sources each bp came from
        var aggregateDataSum = 0.0
        for (bp in 0 until aggregateChrData.size) {
            if (aggregateSources[bp] == 0) continue
            aggregateChrData[bp] /= aggregateSources[bp].toDouble()
            aggregateDataSum += aggregateChrData[bp]
        }

        log.info { "Peaks-only aggregate Created! Smoothing aggregate pile-up for $chr..." }
        return PileUp(aggregateChrData, chr, length, aggregateDataSum)
    }

    override fun peaks(pdf: PDF): List<Region> {
        val allChrPeaks = checkNotNull(allCurrentChrPeaks) { "Top-Down peaks called before pile-up" }
            .map { it.subList(0,20) }
        return mergePeaks(allChrPeaks)
    }

}