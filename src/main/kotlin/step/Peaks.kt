package step

import com.google.common.collect.Iterators
import com.google.common.collect.PeekingIterator
import model.Region
import mu.KotlinLogging

private val log = KotlinLogging.logger {}

/**
 * Find the peaks for the previous calculated PDF over a given threshold.
 *
 * @param threshold How many standard deviations above average the pdf value needs to be to serve
 * as a cut-off for peaks
 */
fun callPeaks(pdf: PDF, threshold: Double): List<Region> {
    val peaks = mutableListOf<Region>()
    var currentRegionStart: Int? = null
    for (chrIndex in 0 until pdf.chrLength) {
        val value = pdf[chrIndex]
        val stdDevsValue = (value - pdf.background.average) / pdf.background.stdDev
        val aboveThreshold =  stdDevsValue > threshold

        if (aboveThreshold && currentRegionStart == null) {
            currentRegionStart = chrIndex
        }
        if(!aboveThreshold && currentRegionStart != null) {
            peaks += Region(currentRegionStart, chrIndex-1)
            currentRegionStart = null
        }
    }

    return peaks
}

fun mergePeaks(vararg allPeaks: List<Region>) = mergePeaks(allPeaks.toList())

/**
 * Merge many lists of peaks into one list of peaks. Peaks that overlap will be merged together.
 */
fun mergePeaks(allPeaks: List<List<Region>>): List<Region> {
    val allPeaksIterators = allPeaks
        .map { peaks -> Iterators.peekingIterator(peaks.sortedBy { it.start }.iterator()) }
    val mergedPeaks = mutableListOf<Region>()

    while (allPeaksIterators.any { it.hasNext() }) {
        var currentRegion = popLowestRegion(allPeaksIterators)
        do {
            var allAboveCurrentRegion = true
            for (peaksIter in allPeaksIterators) {
                if (peaksIter.hasNext() && peaksIter.peek().start <= currentRegion.end) {
                    val nextPeak = peaksIter.next()
                    if (nextPeak.end > currentRegion.end) currentRegion = currentRegion.copy(end = nextPeak.end)
                    allAboveCurrentRegion = false
                }
            }
        } while (!allAboveCurrentRegion)
        mergedPeaks += currentRegion
    }
    return mergedPeaks
}

private fun popLowestRegion(peaksIterators: List<PeekingIterator<Region>>): Region {
    var lowestIter: PeekingIterator<Region>? = null
    for (peaksIter in peaksIterators) {
        if (peaksIter.hasNext() && (lowestIter == null || peaksIter.peek().start < lowestIter.peek().start)) {
            lowestIter = peaksIter
        }
    }
    return lowestIter!!.next()
}