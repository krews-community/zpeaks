package runner

import step.PileUp
import step.runPileUp

class SingleFileZRunner(runConfig: ZRunConfig) : ZRunner("Single File", runConfig) {

    override fun pileUp(chr: String, length: Int, onRange: IntRange?, subsetSize: Int?): PileUp = with(runConfig) {
        val pileUpInput = pileUpInputs[0]
        return runPileUp(pileUpInput.bam, chr, length, pileUpInput.options)
    }
}