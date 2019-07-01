package runner

import model.*
import step.*
import mu.KotlinLogging

private val log = KotlinLogging.logger {}

class TopDownZRunner(runConfig: ZRunConfig) : ZRunner("Top-Down", runConfig) {

    private var allCurrentChrPeaks: List<List<Region>>? = null

    override fun pileUp(chr: String, chrLength: Int, range: IntRange): PileUp = with(runConfig) {
        log.info { "Creating peaks-only aggregate track..." }
        // Array containing all aggregated data
        val aggregateChrData = FloatArray(chrLength)
        // Array containing the number of sources the data came from at each BP
        val aggregateSources = IntArray(chrLength)
        val allPeaks = mutableListOf<List<Region>>()
        allCurrentChrPeaks = allPeaks

        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, chrLength, range, pileUpInput.options)
            val pdf = pdf(pileUp, smoothing)
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
        for (i in 0 until aggregateChrData.size) {
            if (aggregateSources[i] == 0) continue
            aggregateChrData[i] /= aggregateSources[i].toFloat()
            aggregateDataSum += aggregateChrData[i]
        }

        log.info { "Peaks-only aggregate Created! Smoothing aggregate pile-up for $chr..." }
        return PileUp(aggregateChrData, chr, chrLength, range, aggregateDataSum)
    }

    override fun peaks(pdf: PDF): List<Region> {
        val allChrPeaks = checkNotNull(allCurrentChrPeaks) { "Top-Down peaks called before pile-up" }
        return mergePeaks(allChrPeaks)
    }

}