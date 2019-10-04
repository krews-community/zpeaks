package step

import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import model.ChromBounds
import model.PileUpInput
import model.Strand
import mu.KotlinLogging
import java.nio.file.Path


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

/**
 * @param strand: Strand that we want to count using pile-up algorithm
 * @param pileUpAlgorithm: Algorithm we use to choose the values that we pile up.
 * @param forwardShift: Shifts forward (plus) strand by this amount.
 * @param reverseShift: Shifts reverse (minus) strand by this amount.
 * @param filter: Only pile up reads which evaluate to true when passed to this function.
 */
data class PileUpOptions(
    val strand: Strand,
    val pileUpAlgorithm: PileUpAlgorithm,
    val forwardShift: Int = 0,
    val reverseShift: Int = 0,
    val filter: (record: SAMRecord) -> Boolean = { _ -> true }
)

data class BamPileUpRunner(
    val pileUpInputs: List<PileUpInput>
) : PileUpRunner {

    private var prepped = false
    private fun ensurePrepped() {
        if (!prepped) {
            getChromsWithBounds()
        }
    }

    override fun getChromsWithBounds(chrFilter: Map<String, IntRange?>?): Map<String, ChromBounds> {
        val ret = prepBams(pileUpInputs.map { it.bam }, chrFilter)
        prepped = true
        return ret
    }

    override fun getPileUps(chr: String, chrLength: Int, range: IntRange): Sequence<PileUp> = sequence {
        ensurePrepped()
        for (pileUpInput in pileUpInputs) {
            yield(runPileUp(pileUpInput.bam, chr, chrLength, range, pileUpInput.options))
        }
    }

}

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
            if (record.referenceName == "*") return@forEach
            if (!options.filter(record)) return@forEach

            var start = pileUpStart(record, pileUpBounds, options.forwardShift, options.reverseShift)
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
                    val readEnd = pileUpEnd(record, pileUpBounds, options.forwardShift, options.reverseShift)
		    val end = if (readEnd > start) readEnd else start
		    if (readEnd < start) start = readEnd
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
        this < bounds.first -> bounds.first
        this > bounds.last -> bounds.last
        else -> this
    }
