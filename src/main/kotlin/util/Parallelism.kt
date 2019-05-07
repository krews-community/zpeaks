package util

import mu.KotlinLogging
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger

private val log = KotlinLogging.logger {}

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

fun <T, R> runParallel (name: String,
                        values: List<T>,
                        parallelism: Int,
                        logProgress: Boolean = true,
                        runForValue: (value: T) -> R): Set<R> {
    log.info { "$name runs initiated with parallelism = $parallelism" }
    val executor = Executors.newFixedThreadPool(parallelism)
    val logExecutor = Executors.newSingleThreadScheduledExecutor()

    val runsInProgress = AtomicInteger()
    val runsFailed = AtomicInteger()
    val runsSucceeded = AtomicInteger()

    // Setup periodic job to log progress
    val logJob = if (logProgress) {
        logExecutor.scheduleAtFixedRate({
            val inProgress = runsInProgress.get()
            val failed = runsFailed.get()
            val succeeded = runsSucceeded.get()
            val remaining = values.size - succeeded - failed
            log.info { "$name in progress: $inProgress, failed: $failed, Succeeded: $succeeded, remaining: $remaining" }
        }, 5, 5, TimeUnit.SECONDS)
    } else null

    val results = mutableSetOf<R>()
    fun writeResult(value: R) = synchronized(results) { results.add(value) }

    // Run each search result into a "Callable" thread object to write the results.
    val tasks = values.map { value ->
        Callable<Unit> {
            runsInProgress.incrementAndGet()
            var failed = false
            try {
                val result = runForValue(value)
                writeResult(result)
            } catch (e: Exception) {
                log.error(e) { "Error occurred during $name for $value" }
                failed = true
            }
            if (failed) runsFailed.incrementAndGet() else runsSucceeded.incrementAndGet()
            runsInProgress.decrementAndGet()
        }
    }

    // Run all the tasks and block until complete
    executor.invokeAll(tasks)

    // Cancel progress logging job
    logJob?.cancel(false)

    // Shutdown ExecutorServices to prevent application termination issues
    executor.shutdown()
    logExecutor.shutdown()

    log.info { "$name runs complete!" }
    return results
}