import model.*
import mu.KotlinLogging
import step.*
import step.subpeaks.fitSkew
import java.nio.file.Path

private val log = KotlinLogging.logger {}

fun run(
    sam: Path, smoothing: Double, threshold: Double, strand: Strand, atacMode: Boolean,
    pileUpAlgorithm: PileUpAlgorithm
) {
    log.info { "Performing pile-up for sam file $sam with strand=$strand, atacMode=$atacMode, pileUpAlgorithm=$pileUpAlgorithm" }
    val pileUps = pileUpSam(sam, strand, atacMode, pileUpAlgorithm)

    val pileUpSummary = pileUps.entries.joinToString("\n") {
        "chromosome ${it.key}: len=${it.value.chromosomeLength} sum=${it.value.sum}"
    }
    log.info { "Pile-up complete with results: \n$pileUpSummary" }

    for ((chr, pileUp) in pileUps) {
        log.info { "Calculating PDF for chromosome $chr pileup data..." }
        val pdf = pdf(chr, pileUp.values, pileUp.sum, pileUp.chromosomeLength, smoothing, true)
        log.info { "Chromosome $chr PDF completed with background ${pdf.background}" }

        if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
            log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
            continue
        }

        val peaks = callPeaks(pdf, threshold)
        log.info { "Peaks called for threshold $threshold. Peaks found: ${peaks.size}" }
        for (peak in peaks) {
            val region = peak.region
            val peakValues = (region.start..region.end).map { x -> pdf.values.getValue(x) }
            val subPeaks = fitSkew(peakValues, region.start)
        }
    }
}
