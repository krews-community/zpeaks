package util

import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

private val log = KotlinLogging.logger {}

const val SQRT2PI = 2.506628275

fun logProgress(name: String, total: Int, run: (tracker: AtomicInteger) -> Unit) {
    val progressTracker = AtomicInteger(0)
    val logExecutor = Executors.newSingleThreadScheduledExecutor()
    val logJob = logExecutor.scheduleAtFixedRate({
        val complete = progressTracker.get()
        val percentage = "%.2f".format(complete.toDouble() / total * 100)
        log.info { "$name - $percentage% complete." }
    }, 5, 5, TimeUnit.SECONDS)
    run(progressTracker)
    logJob?.cancel(false)
    logExecutor.shutdown()
}

/**
 * Convenience method to util.increment for a key on a mutable map of integers.
 */
fun <T> MutableMap<T, Int>.increment(key: T, incrementBy: Int = 1, defaultValue: Int = 0) {
    this.putIfAbsent(key, defaultValue)
    this[key] = this.getValue(key) + incrementBy
}

/**
 * Convenience method to util.increment for a key on a mutable map of doubles.
 */
fun <T> MutableMap<T, Double>.increment(key: T, incrementBy: Double = 1.0, defaultValue: Double = 0.0) {
    this.putIfAbsent(key, defaultValue)
    this[key] = this.getValue(key) + incrementBy
}

/**
 * @returns the index of the minimum value within the given range
 */
fun <T : Comparable<T>> List<T>.indexOfMin(withinRange: IntRange? = null): Int? {
    val withinList = if (withinRange != null) this.subList(withinRange.first, withinRange.last) else this
    return withinList.withIndex().minBy { it.value }?.index?.plus(withinRange?.first ?: 0)
}

fun splitAtMin(values: List<Double>): Pair<List<Double>, List<Double>> {
    // look only within the middle half to make sure we shrink the sizes by a reasonable amount
    val quarterSize = values.size / 4
    val minIndex = values.subList(quarterSize, quarterSize * 3).indexOfMin()!! + quarterSize

    // Return two halves, split at min index.
    return values.subList(0, minIndex) to values.subList(minIndex, values.size)
}

fun Double.pow(x: Int) = FastMath.pow(this, x)
fun Int.pow(x: Int) = FastMath.pow(this.toDouble(), x)