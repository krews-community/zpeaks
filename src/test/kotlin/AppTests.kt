import model.*
import mu.KotlinLogging
import org.junit.jupiter.api.*
import step.*
import step.subpeaks.scaleSpaceZeros
import util.*

private val log = KotlinLogging.logger {}

class AppTests {

    @Disabled @Test
    fun `Test App Run`() {
        val testBamPath = getResourcePath("ENCFF375IJW.chr22.bam")
        run(testBamPath, 50.0, 6.0, Strand.BOTH, false, PileUpAlgorithm.START)
    }

    @Test
    fun `test candidate selection`() {
        val dist = gaussianDistribution(5.0, 10.0, 15.0, 50)
        val firstHalf = (0 until 25).map { dist[25-it] }
        val secondHalf = (25 until 50).map { dist[it-25] }
        val values = firstHalf + secondHalf
        val candidates = scaleSpaceZeros(values.toDoubleArray())
        log.info { "values: $values" }
        log.info { "candidates: $candidates" }
    }

}