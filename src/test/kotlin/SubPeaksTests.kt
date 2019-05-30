import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.*
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import step.subpeaks.*

private val log = KotlinLogging.logger {}

class SubPeaksTests {

    @Test
    fun `test candidate region selection`() {
        val values = sampleGaussian(4.0, 2.0, 25, 50)

        val candidateRegions = findCandidates(values)
        assertThat(candidateRegions.size).isEqualTo(1)

        val candidateGaussians = candidateGaussians(values, candidateRegions, SkewFitter::initParameters)
        assertThat(candidateGaussians.size).isEqualTo(1)
        assertThat(candidateGaussians[0].parameters.mean)
    }
}

private fun sampleGaussian(u: Double, a: Double, x: Int, size: Int): List<Double> {
    return (0 until size).map { a * exp( -pow( (it - x) / u, 2.0) / 2.0) }
}