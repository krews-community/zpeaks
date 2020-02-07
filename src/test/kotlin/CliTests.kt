import io.SignalOutputFormat
import io.mockk.*
import model.*
import org.junit.jupiter.api.Test
import runner.*
import step.*
import util.*
import java.nio.file.Files

class CliTests {

    @Test
    fun `Test CLI`() {
        val testDir = Files.createTempDirectory("zpeaks_test")
        val signalOut = testDir.resolve("testSignalOut.bed")
        val peaksOut = testDir.resolve("testPeaksOut.bed")

        val mockRun = spyk<(RunType, ZRunConfig, Boolean) -> Unit>()
        val args =
            """
            -bamIn=$MULTI_BAM_1_PATH -strand=both -pileUpAlgorithm=start -forwardShift=5 -reverseShift=10
            -bamIn=$MULTI_BAM_2_PATH -strand=plus -pileUpAlgorithm=length -forwardShift=-5 -reverseShift=15
            -chrFilter=$CHR_FILTER_PATH
            -signalOut=$signalOut -signalOutType=raw -signalOutFormat=wig -signalResolution=2
            -peaksOut=$peaksOut
            -smoothingFactor=60.0 -threshold=7.0
            -parallelism=4
            """.trimIndent().split("\\s+".toRegex())
        ZPeaksCommand(run = mockRun).main(args)

        fun matchPileUpInput(it: PileUpInput, other: PileUpInput): Boolean {	
	    return it.bam.equals(other.bam) && it.options.strand.equals(other.options.strand) &&
                it.options.pileUpAlgorithm.equals(other.options.pileUpAlgorithm) &&
		it.options.forwardShift == other.options.forwardShift &&
		it.options.reverseShift == other.options.reverseShift
        }
	
        fun matchBamPileUpRunner(it: BamPileUpRunner, other: BamPileUpRunner): Boolean {
	    if (it.pileUpInputs.size != other.pileUpInputs.size) return false
	    var matched = true	
	    other.pileUpInputs.forEachIndexed { i, v ->
	        if (!matchPileUpInput(it.pileUpInputs[i], v)) matched = false
	    }
	    return matched
        }

        verify { mockRun(any(), match {
	    it.chrFilter != null && it.chrFilter!!.equals(mapOf("chr22" to (30000000 until 30500000))) &&
	    it.peaksOut != null && it.peaksOut!!.equals(peaksOut) &&
	    it.smoothing == 60.0 && it.threshold == 7.0 && it.parallelism == 4 &&
  	    it.signalOut != null && it.signalOut!!.equals(
	        SignalOutput(signalOut, SignalOutputType.RAW, SignalOutputFormat.WIG, signalResolution = 2)
	    ) &&
	    it.pileUpRunner is BamPileUpRunner && matchBamPileUpRunner(
		it.pileUpRunner as BamPileUpRunner, BamPileUpRunner(listOf(
                    PileUpInput(MULTI_BAM_1_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START, 5, 10)),
                    PileUpInput(MULTI_BAM_2_PATH, PileUpOptions(Strand.PLUS, PileUpAlgorithm.LENGTH, -5, 15))
                ))
	    )
	}, false) }
    }

}
