package runner

import mu.KotlinLogging
import step.*

private val log = KotlinLogging.logger {}

class BottomUpZRunner(runConfig: ZRunConfig) : ZRunner("Bottom-Up", runConfig) {

    override fun pileUp(chr: String, length: Int, onRange: IntRange?, subsetSize: Int?): PileUp = with(runConfig) {
        log.info { "Creating aggregate pile-up for $chr..." }
        val aggregatePileUpData = DoubleArray(length)
        var sum = 0.0
        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, length, pileUpInput.options)
            for (index in 0 until length) {
                val value = pileUp.scaledValue(index)
                aggregatePileUpData[index] += value
                sum += value
            }
        }

        log.info { "Aggregate pile-up for $chr complete!" }
        return PileUp(aggregatePileUpData, chr, length, sum)
    }

}