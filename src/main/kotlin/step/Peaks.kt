package step

import model.Region
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.*

private val log = KotlinLogging.logger {}

data class Peak(val region: Region, val score: Double)

fun callPeaks(pdfs: Map<String, PDF>, threshold: Double): Map<String, List<Peak>> {
    log.info { "Calling peaks..." }
    val peaks = mutableMapOf<String, List<Peak>>()
    for ((chr, pdf) in pdfs) {
        peaks[chr] = callChromPeaks(pdf, threshold)
    }
    log.info { "Peak calling complete!" }
    return peaks
}

/**
 * Find the peaks for the previous calculated PDF over a given threshold.
 *
 * @param threshold How many standard deviations above average the pdf value needs to be to serve
 * as a cut-off for peaks
 */
fun callChromPeaks(pdf: PDF, threshold: Double): List<Peak> {
    val peaks = mutableListOf<Peak>()
    var currentRegionStart: Int? = null
    var currentRegionMax = 0.0
    for (chrIndex in 0 until pdf.chrLength) {
        val value = pdf[chrIndex]
        currentRegionMax = max(currentRegionMax, value)
        val stdDevsValue = (value - pdf.background.average) / pdf.background.stdDev
        val aboveThreshold =  stdDevsValue > threshold

        if (aboveThreshold && currentRegionStart == null) {
            currentRegionStart = chrIndex
        }
        if(!aboveThreshold && currentRegionStart != null) {
            peaks += Peak(Region(currentRegionStart, chrIndex-1), currentRegionMax)
            currentRegionStart = null
            currentRegionMax = 0.0
        }
    }

    return peaks
}