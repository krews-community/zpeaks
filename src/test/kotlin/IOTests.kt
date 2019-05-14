import model.Strand
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import step.PileUpAlgorithm
import util.*
import java.nio.file.*

class IOTests {

    private fun testPileUpFormat(fileName: String, format: PileUpOutputFormat, checkContents: Boolean = true) {
        val outputPath = Files.createTempDirectory("zpeaks_test").resolve(fileName)
        runPileUp(TEST_BAM_PATH, Strand.BOTH, false, PileUpAlgorithm.LENGTH,
            PileUpOutput(outputPath, format))
        //Files.copy(outputPath, testBamPath.resolveSibling(fileName))
        if (checkContents) assertThat(outputPath).hasSameContentAs(getResourcePath(fileName))
    }

    @Test
    fun `test pile-up output in wig format`() {
        testPileUpFormat("ENCFF375IJW.chr22.pileup.wig", PileUpOutputFormat.WIG)
    }

    @Test
    fun `test pile-up output in bigwig format`() {
        testPileUpFormat("ENCFF375IJW.chr22.pileup.bw", PileUpOutputFormat.BIG_WIG, checkContents = false)
    }

    @Test
    fun `test pile-up output in bed-graph format`() {
        testPileUpFormat("ENCFF375IJW.chr22.pileup.bedGraph", PileUpOutputFormat.BED_GRAPH)
    }

}