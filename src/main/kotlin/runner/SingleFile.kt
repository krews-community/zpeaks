package runner

import step.*

class SingleFileZRunner(runConfig: ZRunConfig) : ZRunner("Single File", runConfig) {

    override fun pileUp(chr: String, chrLength: Int, range: IntRange): PileUp = with(runConfig) {
        return pileUpRunner.getPileUps(chr, chrLength, range).toList()[0]
    }
}