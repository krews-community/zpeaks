FROM openjdk:12-jdk-alpine as build
COPY . /src
RUN ./gradlew clead shadowJar

FROM openjdk:12-jre-alpine
COPY --from=build /src/build/*.jar /app/zpeaks.jar