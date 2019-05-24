package util

import com.google.common.base.CaseFormat
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath
import java.util.concurrent.*
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

fun <T> runParallel (name: String, values: List<T>, runForValue: (value: T) -> Unit) {
    val commonPool = ForkJoinPool.commonPool()
    log.info { "$name initiated: ${values.size} runs with parallelism = ${commonPool.parallelism}" }

    logProgress(name, values.size) { tracker ->
        // Run each search result into a "Callable" thread object to write the results.
        val tasks = values.map { value ->
            Callable<Unit> {
                try {
                    runForValue(value)
                } catch (e: Exception) {
                    log.error(e) { "Error during $name" }
                }
                tracker.incrementAndGet()
            }
        }

        // Run all the tasks and block until complete
        commonPool.invokeAll(tasks)
    }

    // Shutdown ExecutorService to prevent application termination issues
    log.info { "$name runs complete!" }
}

/**
 * Extension function to get enum name in lower-hyphen-format
 */
val Enum<*>.lowerHyphenName: String get() = CaseFormat.UPPER_UNDERSCORE.to(CaseFormat.LOWER_HYPHEN, this.name)

fun Double.pow(x: Int) = FastMath.pow(this, x)