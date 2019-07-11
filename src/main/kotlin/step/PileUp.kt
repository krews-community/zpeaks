package step

import htsjdk.samtools.*
import model.*
import mu.KotlinLogging
import util.*
import java.nio.file.Path
import kotlin.math.sqrt


private val log = KotlinLogging.logger {}

enum class PileUpAlgorithm {
    // Use only the start of each alignment
    START,
    // Use the mid-point between start and end for each alignment. *Paired End Only
    MID_POINT,
    // Use every BP between start and end for each alignment. *Paired End Only
    LENGTH
}

// When using the LENGTH pile-up algorithm we will not include values for lengths that are larger than this value.
// They are most likely junk and will waste resources.
const val LENGTH_LIMIT = 100_000

class PileUp(
    // Pile-up values in chromosome.
    private val values: FloatArray,
    // Chromosome as named in alignment file
    override val chr: String,
    // Chromosome length pulled directly from BAM File
    val chrLength: Int,
    // Chromosome range to be used downstream
    override val range: IntRange,
    // Sum of all pile-up values in chromosome calculated on the fly and cached for efficiency
    val sum: Double
) : SignalData {
    override operator fun get(bp: Int): Float = values[bp]

    private val scalingFactor: Float = if(sum == 0.0) 1.0F else sqrt(sum).toFloat()
    fun scaledValue(bp: Int): Float = this[bp] / scalingFactor
}

/**
 * @param strand: Strand that we want to count using pile-up algorithm
 * @param pileUpAlgorithm: Algorithm we use to choose the values that we pile up.
 * @param forwardShift: shifts forward (plus) strand by this amount
 * @param reverseShift: shifts reverse (minus) strand by this amount
 */
data class PileUpOptions(
    val strand: Strand,
    val pileUpAlgorithm: PileUpAlgorithm,
    val forwardShift: Int = 0,
    val reverseShift: Int = 0
)

/**
 * Reads a SAM or BAM file and creates a pile-up representation in memory.
 *
 * @return the in-memory pile-up
 */
fun runPileUp(bam: Path, chr: String, chrLength: Int, range: IntRange, options: PileUpOptions): PileUp {
    log.info { "Performing pile-up for chromosome $chr on alignment $bam..." }
    val values = FloatArray(chrLength) { 0.0F }
    var sum = 0L
    val pileUpBounds = 0 until chrLength
    SamReaderFactory.make().validationStringency(ValidationStringency.SILENT).open(bam).use { reader ->
        reader.query(chr, 0, 0, false).forEach { record ->
            // If we're only using plus strand values and this record is a plus strand and
            // visa versa for minus strand, continue.
            if ((options.strand == Strand.PLUS && record.readNegativeStrandFlag) ||
                (options.strand == Strand.MINUS && !record.readNegativeStrandFlag)
            ) return@forEach

            // If the record is junk, continue
            if (record.end < record.start || record.end - record.start > LENGTH_LIMIT) return@forEach
	    if (record.referenceName.equals("*")) return@forEach

            val start = pileUpStart(record, pileUpBounds, options.forwardShift, options.reverseShift)
            when (options.pileUpAlgorithm) {
                PileUpAlgorithm.START -> {
                    values[start]++
                    sum++
                }
                PileUpAlgorithm.MID_POINT -> {
                    val end = pileUpEnd(record, pileUpBounds, options.forwardShift, options.reverseShift)
                    val midPoint = (start + end) / 2
                    values[midPoint]++
                    sum++
                }
                PileUpAlgorithm.LENGTH -> {
                    val end = pileUpEnd(record, pileUpBounds, options.forwardShift, options.reverseShift)
                    val length = end - start
                    for (i in start until end) {
                        values[i]++
                    }
                    sum += length
                }
            }
        }
    }

    log.info { "Pile-up complete with sum $sum" }
    return PileUp(values, chr, chrLength, range, sum.toDouble())
}

/**
 * Calculate the "start" value for our pile-up algorithm for the given record and settings.
 *
 * This start will actually be the record's end it's strand is minus
 */
private fun pileUpStart(record: SAMRecord, bounds: IntRange, forwardShift: Int, reverseShift: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        (record.start + forwardShift).withinBounds(bounds)
    } else {
        (record.end + reverseShift).withinBounds(bounds)
    }
}

/**
 * Calculate the "end" value for our pile-up algorithm for the given SAM record and settings.
 *
 * This end will actually be the record's start it's strand is minus
 */
private fun pileUpEnd(record: SAMRecord, bounds: IntRange, forwardShift: Int, reverseShift: Int): Int {
    return if (!record.readNegativeStrandFlag) {
        (record.end + reverseShift).withinBounds(bounds)
    } else {
        (record.start + forwardShift).withinBounds(bounds)
    }
}

/**
 * @returns: if < start: start, if > end: end, else value
 */
private fun Int.withinBounds(bounds: IntRange) =
    when {
        this < bounds.start -> bounds.start
        this > bounds.endInclusive -> bounds.endInclusive
        else -> this
    }