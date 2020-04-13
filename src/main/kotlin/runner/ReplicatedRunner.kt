package runner

import io.*
import model.FitMode
import model.Region
import model.SignalOutputType
import model.ReplicatedRegion
import mu.KotlinLogging
import step.*
import step.subpeaks.SkewFitter
import step.subpeaks.SkewReplicatedFitter
import step.subpeaks.StandardFitter
import step.subpeaks.StandardReplicatedFitter
import java.util.concurrent.ForkJoinPool

private val log = KotlinLogging.logger {}

class ReplicatedRunner(private val runConfig: ZRunConfig) {

    fun pileUp(chr: String, chrLength: Int, range: IntRange): List<PileUp> = with(runConfig) {
        log.info { "Creating aggregate pile-up for $chr..." }
        val pileUps: MutableList<PileUp> = mutableListOf()
        pileUpRunner.getPileUps(chr, chrLength, range).forEach { pileUps.add(it) }
        return pileUps
    }

    fun pdf(pileUp: List<PileUp>): List<PDF> = with(runConfig) { pileUp.map { pdf(it, smoothing) } }

    fun peaks(pdf: List<PDF>): List<Region> = with(runConfig) {
        val unmergedPeaks: List<ReplicatedRegion> = pdf.map { callPeaks(it, threshold) }.foldIndexed(listOf(), { i, acc, it ->
            acc + it.map { r -> ReplicatedRegion(r.start, r.end, setOf(i)) }
        })
        val sortedPeaks = unmergedPeaks.sortedBy { it.start }
	if (sortedPeaks.size === 0) return mutableListOf()
        var current = sortedPeaks[0]
        val ret: MutableList<Region> = mutableListOf()
        sortedPeaks.forEach {
            current = if (current.end > it.start)
                ReplicatedRegion(
		    start = current.start,
		    end = if (current.end < it.end) it.end else current.end,
		    replicates = current.replicates union it.replicates
		)
            else {
                if (current.replicates.size === pdf.size) ret.add(Region(current.start, current.end))
                it
            }
        }
        return ret
    }

    fun run() = with(runConfig) {

        if (parallelism != null)
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", parallelism.toString())
        log.info { "ZPeaks run started with parallelism = ${ForkJoinPool.commonPool().parallelism}" }
        log.info { "Running replicated ZPeaks operation" }

        val chrsWithBounds = pileUpRunner.getChromsWithBounds(chrFilter)
        for ((chr, bounds) in chrsWithBounds) {

            val pileUp = pileUp(chr, bounds.length, bounds.range)
            val pdf = pdf(pileUp)
            var zero = false
            pdf.forEach {
                if (0.0 == it.background.average || 0.0 == it.background.stdDev) {
                    log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
                    zero = true
                }
            }
            if (zero) continue

            if (signalOut != null)
                log.warn {
                    "signal file generation is not supported for replicated runs; run ZPeaks separately for each file"
                }
            if (peaksOut == null) continue

            val peaks = peaks(pdf)
            if (fitMode == FitMode.SKEW) {
                val skewSubPeaks = SkewReplicatedFitter.fit(chr, peaks, pdf)
                writeReplicatedSkewSubPeaksBed(peaksOut, chr, skewSubPeaks)
            } else {
                val subPeaks = StandardReplicatedFitter.fit(chr, peaks, pdf)
                writeReplicatedStandardSubPeaksBed(peaksOut, chr, subPeaks)
            }

        }
        log.info { "Replicated ZPeaks run complete!" }
    }

}
