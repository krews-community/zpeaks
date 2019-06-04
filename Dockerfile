FROM openjdk:12-jdk-alpine as build
COPY . /src
WORKDIR /src
RUN /src/gradlew clean shadowJar

FROM openjdk:12-jre-alpine
ENV JAVA_OPTS="-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap"
COPY --from=build /src/build/*.jar /app/zpeaks.jar