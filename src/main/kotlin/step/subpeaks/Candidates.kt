package step.subpeaks

import model.Region
import org.apache.commons.math3.util.FastMath
import util.scaleSpaceSmooth

/*
 * Uses the scale-space image of the data to generate an initial list
 * of candidate Gaussians. Zero crossings as they emerge from most blurred
 * to least define the Gaussians' parameters.
 *
 * values: input vector from which to generate the scale space image.
 */
fun findCandidates(values: List<Double>): List<Region> {
    var candidates = mutableListOf<Region>()

    // Start with a kernel near the square root of the vector size and loop down to 2
    val kernel = Math.floor(FastMath.sqrt(values.size.toDouble()) * 1.1).toInt()
    for (kernelSize in kernel downTo 2 step 1) {

        // Smooth at this kernel width
        val blurred = scaleSpaceSmooth(values, kernelSize.toDouble())

        // Find zero crossings
        val zeroCrossings = mutableListOf<Int>()
        var cSign = blurred[0] > 0
        for (j in 1 until blurred.size) {
            if (cSign != blurred[j] > 0) {
                cSign = blurred[j] > 0
                zeroCrossings += j
            }
        }

        // Find matching zero-crossings in existing candidates
        val currentCandidates = mutableListOf<Region>()
        val matchedZeroCrossings = mutableSetOf<Int>()
        for (candidate in candidates) {
            val candidateLength = candidate.end - candidate.start
            // Find the closest matching zero crossings within candidateLength to the current candidate
            var startMatch: Int? = null
            var startMatchDist: Int? = null
            var endMatch: Int? = null
            var endMatchDist: Int? = null
            for (crossing in zeroCrossings) {
                if (matchedZeroCrossings.contains(crossing)) continue
                val startDist = FastMath.abs(crossing - candidate.start)
                val endDist = FastMath.abs(crossing - candidate.end)
                if (startDist < candidateLength && (startMatch == null || startDist < startMatchDist!!)) {
                    startMatch = crossing
                    startMatchDist = startDist
                } else if (endDist < candidateLength && (endMatch == null || endDist < endMatchDist!!)) {
                    endMatch = crossing
                    endMatchDist = endDist
                }
            }

            // If we got a match, add new match to candidates in place of old. Otherwise add old candidate.
            if (startMatch != null && endMatch != null) {
                matchedZeroCrossings += startMatch
                matchedZeroCrossings += endMatch
                currentCandidates += Region(startMatch, endMatch)
            } else {
                currentCandidates += candidate
            }
        }

        // Add unmatched zero crossings as new candidates
        val unmatchedZeroCrossings = zeroCrossings.filter { !matchedZeroCrossings.contains(it) }
        for (n in 1 until unmatchedZeroCrossings.size step 2) {
            currentCandidates += Region(unmatchedZeroCrossings[n - 1], unmatchedZeroCrossings[n])
        }

        candidates = currentCandidates
    }

    /*
    val filtered = candidates.filter { c1 ->
        var containsAnother = false
        for(c2 in candidates) {
            if (c1.start < c2.start && c1.end > c2.end) {
                log.info { "eliminating $c1. contains $c2" }
                containsAnother = true
                break
            }
        }
        !containsAnother
    }
    */

    return candidates.sortedBy { it.start }
}