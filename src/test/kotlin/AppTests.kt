import model.*
import org.junit.jupiter.api.*
import step.*
import util.getResourcePath

class AppTests {

    @Disabled @Test
    fun `Test App Run`() {
        val testBamPath = getResourcePath("ENCFF375IJW.chr22.bam")
        run(testBamPath, 50.0, 6.0, Strand.BOTH, false, PileUpAlgorithm.START)
    }

}