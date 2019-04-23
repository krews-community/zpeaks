package step

import java.lang.Math.*
import java.util.*
import model.*
import mu.KotlinLogging
import util.*

private val log = KotlinLogging.logger {}

const val BACKGROUND_LIMIT = 1000

data class Background(val average: Double, val stdDev: Double)
data class PDF(val values: SortedMap<Int, Double>, val background: Background)

fun pdf(pileUpValues: Map<Int, Int>, pileUpSum: Int, chrLength: Int, bandwidth: Double, noAmplitude: Boolean): PDF {
    val windowSize = windowSize(bandwidth)
    val lookupTable = lookupTable(noAmplitude, windowSize, bandwidth, pileUpSum)
    val pdfValues = sortedMapOf<Int, Double>()

    for ((chrIndex, pileUpValue) in pileUpValues) {
        pdfValues.increment(chrIndex, pileUpValue * lookupTable[0])
        for (i in 1 until windowSize) {
            if (chrIndex + i < chrLength) {
                pdfValues.increment(chrIndex + i, pileUpValue * lookupTable[i])
            }
            if (chrIndex - i > 0) {
                pdfValues.increment(chrIndex - i, pileUpValue * lookupTable[i])
            }
        }
    }

    val background = background(lookupTable, pileUpSum, windowSize, chrLength)
    return PDF(pdfValues, background)
}

fun callPeaks(pdf: PDF, threshold: Double): List<Region> {
    val peaks = mutableListOf<Region>()
    var currentRegionStart: Int? = null
    var currentRegionEnd: Int? = null
    var previousChrIndex: Int? = null
    for ((chrIndex, value) in pdf.values) {
        // Since the step.PDF values are stored as a sparse ordered map where any missing values are 0,
        // we need to check if we skipped any values. If so cut off the region and add it to peaks.
        if (previousChrIndex != null && previousChrIndex == chrIndex - 1 && currentRegionStart != null) {
            peaks += Region(currentRegionStart, currentRegionEnd!!)
            currentRegionStart = null
            currentRegionEnd = null
        }

        if (value - pdf.background.average / pdf.background.stdDev > threshold) {
            if (currentRegionStart != null) {
                currentRegionStart = chrIndex
            }
            currentRegionEnd = chrIndex
        } else if(currentRegionStart != null) {
            peaks += Region(currentRegionStart, currentRegionEnd!!)
            currentRegionStart = null
            currentRegionEnd = null
        }

        previousChrIndex = chrIndex
    }

    return peaks
}

private fun windowSize(bandwidth: Double): Int {
    return (sqrt(log(Double.MIN_VALUE * SQRT2PI) * -2.0) * bandwidth).toInt()
}

/**
 * Compute the step.background. We get the average number of regions per window, then get n
 * random distributions of this many reads across the window. We repeat 1000 times; for
 * each repeat, we compute the step.PDF at the center; at high numbers of repeats, the
 * distribution of PDFs should be near normal.
 */
private fun background(lookupTable: List<Double>, numBins: Int, windowSize: Int, chrLength: Int): Background {
    val averageN = numBins * windowSize / chrLength
    val backgroundDist = mutableListOf<Double>()
    for (i in 0 until BACKGROUND_LIMIT) {
        for (j in 0 until averageN) {
            val x = abs(random() % (windowSize + 1) - windowSize / 2).toInt()
            backgroundDist[i] += lookupTable[x]
        }
    }

    // Compute step.background average and standard deviation from the distribution.
    val average = backgroundDist.sum() / BACKGROUND_LIMIT
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