package step

import java.lang.Math.*
import model.*
import util.*
import java.util.concurrent.ConcurrentHashMap
import kotlin.random.Random

const val BACKGROUND_LIMIT = 1000

data class Background(val average: Double, val stdDev: Double)
data class PDF(val values: Map<Int, Double>, val background: Background, val chrLength: Int)

/**
 * Using FSeq algorithm, calculate the probability density function for the given pile-up data
 */
fun pdf(chr: String, pileUpValues: Map<Int, Int>, pileUpSum: Int, chrLength: Int, bandwidth: Double, noAmplitude: Boolean): PDF {
    val windowSize = windowSize(bandwidth)
    val lookupTable = lookupTable(noAmplitude, windowSize, bandwidth, pileUpSum)
    val pdfValues = ConcurrentHashMap<Int, Double>(pileUpValues.size * 2)

    logProgress("FSeq PDF on $chr", pileUpValues.size) { tracker ->
        pileUpValues.entries.parallelStream().forEach { (chrIndex, pileUpValue) ->
            pdfValues.increment(chrIndex, pileUpValue * lookupTable[0])
            for (i in 1 until windowSize) {
                if (chrIndex + i < chrLength) {
                    pdfValues.increment(chrIndex + i, pileUpValue * lookupTable[i])
                }
                if (chrIndex - i > 0) {
                    pdfValues.increment(chrIndex - i, pileUpValue * lookupTable[i])
                }
            }
            tracker.incrementAndGet()
        }
    }

    val background = background(lookupTable, pileUpSum, windowSize, chrLength)
    return PDF(pdfValues, background, chrLength)
}

data class Peak(val region: Region, val score: Double)

/**
 * Find the peaks for the previous calculated PDF over a given threshold.
 *
 * @param threshold How many standard deviations above average the pdf value needs to be to serve
 * as a cut-off for peaks
 */
fun callPeaks(pdf: PDF, threshold: Double): List<Peak> {
    val peaks = mutableListOf<Peak>()
    var currentRegionStart: Int? = null
    var currentRegionMax = 0.0
    for (chrIndex in 0 until pdf.chrLength) {
        val value = pdf.values.getOrDefault(chrIndex, 0.0)
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

private fun windowSize(bandwidth: Double): Int {
    return (sqrt(log(Double.MIN_VALUE * SQRT2PI) * -2.0) * bandwidth).toInt()
}

/**
 * Compute the background. We get the average number of regions per window, then get n
 * random distributions of this many reads across the window. We repeat 1000 times; for
 * each repeat, we compute the PDF at the center; at high numbers of repeats, the
 * distribution of PDFs should be near normal.
 */
private fun background(lookupTable: List<Double>, numBins: Int, windowSize: Int, chrLength: Int): Background {
    val averageN = numBins * windowSize / chrLength
    val backgroundDist = mutableListOf<Double>()
    for (i in 0 until BACKGROUND_LIMIT) {
        var bgVal = 0.0
        for (j in 0 until averageN) {
            val x = Random.nextInt(0, windowSize / 2)
            bgVal += lookupTable[x]
        }
        backgroundDist += bgVal
    }

    // Compute step.background average and standard deviation from the distribution.
    val average = backgroundDist.average()
    val stdDev = sqrt(backgroundDist.sumByDouble { pow(it - average, 2.0) } / BACKGROUND_LIMIT)
    return Background(average, stdDev)
}

private fun lookupTable(noAmplitude: Boolean, windowSize: Int, bandwidth: Double, n: Int): List<Double> {
    val lookupTable = if (windowSize > 0) {
        val a = if (noAmplitude) 1.0 else 1.0 / bandwidth / SQRT2PI
        gaussianDistribution(a, 0.0, bandwidth, windowSize)
    } else {
        listOf(1.0)
    }
    return if (!noAmplitude) lookupTable.map { it / n } else lookupTable
}