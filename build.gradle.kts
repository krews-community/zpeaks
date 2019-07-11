import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar

plugins {
    kotlin("jvm") version "1.3.21"
    id("com.github.johnrengelman.shadow") version "4.0.2"
    id("application")
}

group = "com.genomealmanac"
version = "1.0.2"
val artifactID = "zpeaks"

repositories {
    jcenter()
    mavenCentral()
}

dependencies {
    implementation(kotlin("stdlib-jdk8"))
    compile("com.github.ajalt", "clikt", "2.0.0")
    compile("io.github.microutils","kotlin-logging","1.6.10")
    compile("ch.qos.logback", "logback-classic","1.2.3")
    compile("com.github.samtools", "htsjdk","2.19.0")
    compile("org.slf4j", "log4j-over-slf4j", "1.7.26")
    compile("org.ejml", "ejml-all", "0.38")
    compile("org.apache.commons", "commons-math3", "3.6.1")
    compile("com.google.guava", "guava", "27.1-jre")

    testImplementation("org.junit.jupiter", "junit-jupiter", "5.4.0")
    testImplementation("io.mockk:mockk:1.9.3")
    testCompile("org.assertj", "assertj-core", "3.11.1")
    testCompile("org.knowm.xchart", "xchart", "3.5.4")
}

tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}

tasks.withType<Test> {
    useJUnitPlatform { }
    testLogging {
        events("passed", "skipped", "failed")
    }
}

application {
    mainClassName = "AppKt"
}
val shadowJar: ShadowJar by tasks
shadowJar.apply {
    baseName = artifactID
    classifier = ""
    destinationDir = file("build")
}
