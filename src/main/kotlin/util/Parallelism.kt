package util

import mu.KotlinLogging
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger

private val log = KotlinLogging.logger {}

fun <T> runParallel (name: String, values: List<T>, parallelism: Int, runForValue: (value: T) -> Unit) {
    log.info { "${values.size} $name runs initiated with parallelism = $parallelism" }
    val executor = Executors.newFixedThreadPool(parallelism)

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
        executor.invokeAll(tasks)
    }

    // Shutdown ExecutorService to prevent application termination issues
    executor.shutdown()
    log.info { "$name runs complete!" }
}