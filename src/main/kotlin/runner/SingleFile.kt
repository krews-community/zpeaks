package runner

import step.*

class SingleFileZRunner(runConfig: ZRunConfig) : ZRunner("Single File", runConfig) {

    override fun pileUp(chr: String, chrLength: Int, range: IntRange): PileUp = with(runConfig) {
        val pileUpInput = pileUpInputs[0]
        return runPileUp(pileUpInput.bam, chr, chrLength, range, pileUpInput.options)
    }
}