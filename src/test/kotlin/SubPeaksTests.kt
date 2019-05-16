import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath.*
import org.assertj.core.api.Assertions.assertThat
import org.junit.jupiter.api.Test
import step.subpeaks.*

private val log = KotlinLogging.logger {}

class SubPeaksTests {

    @Test
    fun `test candidate selection`() {
        val values = sampleGaussian(4.0, 2.0, 25, 50)
        val candidates = findCandidates(values)
        assertThat(candidates.size).isEqualTo(1)
    }
}

private fun sampleGaussian(u: Double, a: Double, x: Int, size: Int): List<Double> {
    return (0 until size).map { a * exp( -pow( (it - x) / u, 2.0) / 2.0) }
}