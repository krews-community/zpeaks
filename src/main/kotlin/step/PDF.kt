package step

import com.google.common.util.concurrent.AtomicDoubleArray
import model.*
import mu.KotlinLogging
import util.*
import java.util.stream.IntStream
import kotlin.random.Random

private val log = KotlinLogging.logger {}

const val BACKGROUND_LIMIT = 1000

data class Background(val average: Double, val stdDev: Double)

class PDF(
    private val values: FloatArray,
    val background: Background,
    override val chr: String,
    override val range: IntRange
): SignalData {
    override operator fun get(bp: Int): Float = values[bp - range.start]
}

/**
 * Using FSeq algorithm, calculate the probability density function for the given pile-up data
 */
fun pdf(pileUp: PileUp, bandwidth: Double): PDF {
    val windowSize = windowSize(bandwidth)
    val lookupTable = lookupTable(windowSize, bandwidth)
    val length = pileUp.range.length
    val pdfValues = AtomicDoubleArray(length)

    logProgress("Creating PDF for ${pileUp.chr}", length) { tracker ->
        IntStream.range(pileUp.range.start, pileUp.range.endInclusive).parallel().forEach { chrIndex ->
            tracker.incrementAndGet()
            val pileUpValue = pileUp[chrIndex]
            if (pileUpValue == 0.0F) return@forEach
            val arrayIndex = chrIndex - pileUp.range.start
            pdfValues.addAndGet(arrayIndex, pileUpValue * lookupTable[0])
            for (i in 1 until windowSize) {
                if (chrIndex + i < pileUp.range.endInclusive) {
                    pdfValues.addAndGet(arrayIndex + i, pileUpValue * lookupTable[i])
                }
                if (chrIndex - i > pileUp.range.start) {
                    pdfValues.addAndGet(arrayIndex - i, pileUpValue * lookupTable[i])
                }
            }
        }
    }

    val activeLength = activeLength(pileUp)
    val background = background(lookupTable, pileUp.sum, windowSize, activeLength)
    // Get values as a FloatArray
    val values = FloatArray(length) { pdfValues[it].toFloat() }

    // Normalize PDF values
    val average = values.average().toFloat()
    for (i in 0 until values.size) values[i] /= average
    val normalizedBackground = Background(background.average / average, background.stdDev / average)

    return PDF(values, normalizedBackground, pileUp.chr, pileUp.range)
}

private fun activeLength(pileUp: PileUp): Int {
    var start: Int? = null
    for (i in 0 until pileUp.chrLength) {
        if (pileUp[i] > 0.0) {
            start = i
            break
        }
    }
    if (start == null) return 0
    var end: Int? = null
    for (i in pileUp.chrLength-1 downTo start) {
        if (pileUp[i] > 0.0) {
            end = i
            break
        }
    }
    return end!! - start
}

private fun windowSize(bandwidth: Double): Int {
    return (Math.sqrt(Math.log(Double.MIN_VALUE * SQRT2PI) * -2.0) * bandwidth).toInt()
}

/**
 * Compute the background. We get the average number of regions per window, then get n
 * random distributions of this many reads across the window. We repeat 1000 times; for
 * each repeat, we compute the PDF at the center; at high numbers of repeats, the
 * distribution of PDFs should be near normal.
 */
private fun background(lookupTable: List<Double>, numBins: Double, windowSize: Int, activeLength: Int): Background {
    if (activeLength == 0) return Background(0.0, 0.0)

    val averageN = numBins * windowSize / activeLength.toFloat()
    val backgroundDist = mutableListOf<Double>()

    if (averageN > 1.0) {
        for (i in 0 until BACKGROUND_LIMIT) {
            var bgVal = 0.0
            for (j in 0 until averageN.toInt()) {
                val x = Random.nextInt(0, windowSize / 2)
                bgVal += lookupTable[x]
            }
            backgroundDist += bgVal
        }
    } else {
        for (i in 0 until BACKGROUND_LIMIT) {
            if (Random.nextDouble() > averageN) continue
            val x = Random.nextInt(0, windowSize / 2)
            backgroundDist += lookupTable[x]
        }
    }

    // Compute step.background average and standard deviation from the distribution.
    val average = backgroundDist.average()
    val stdDev = Math.sqrt(backgroundDist.sumByDouble { Math.pow(it - average, 2.0) } / BACKGROUND_LIMIT)
    return Background(average, stdDev)
}

private fun lookupTable(windowSize: Int, bandwidth: Double): List<Double> {
    return if (windowSize > 0) {
        gaussianDistribution(1.0, 0.0, bandwidth, windowSize)
    } else {
        listOf(1.0)
    }
}