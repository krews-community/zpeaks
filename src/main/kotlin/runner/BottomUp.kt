package runner

import mu.KotlinLogging
import step.*

private val log = KotlinLogging.logger {}

class BottomUpZRunner(runConfig: ZRunConfig) : ZRunner("Bottom-Up", runConfig) {

    override fun pileUp(chr: String, chrLength: Int, range: IntRange): PileUp = with(runConfig) {
        log.info { "Creating aggregate pile-up for $chr..." }
        val aggregatePileUpData = FloatArray(chrLength)
        var sum = 0.0
        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, chrLength, range, pileUpInput.options)
            for (index in 0 until chrLength) {
                val value = pileUp.scaledValue(index)
                aggregatePileUpData[index] += value
                sum += value
            }
        }

        log.info { "Aggregate pile-up for $chr complete!" }
        return PileUp(aggregatePileUpData, chr, chrLength, range, sum)
    }

}