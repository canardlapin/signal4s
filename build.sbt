import org.scalajs.linker.interface.ModuleKind
import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import sbtcrossproject.CrossPlugin.autoImport.*
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

ThisBuild / organization := "io.github.canardlapin"
ThisBuild / scalaVersion := "3.7.4"
ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / versionScheme := Some("early-semver")
ThisBuild / homepage := Some(url("https://github.com/canardlapin/signal4s"))
ThisBuild / licenses := Seq(
  "Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt")
)
ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/canardlapin/signal4s"),
    "scm:git:https://github.com/canardlapin/signal4s.git",
    Some("scm:git:git@github.com:canardlapin/signal4s.git")
  )
)
ThisBuild / developers := List(
  Developer(
    id = "canardlapin",
    name = "canardlapin",
    email = "307091466+canardlapin@users.noreply.github.com",
    url = url("https://github.com/canardlapin")
  )
)
ThisBuild / scalacOptions ++= Seq(
  "-deprecation",
  "-feature",
  "-unchecked",
  "-Wunused:all",
  "-Wvalue-discard",
  "-Werror"
)
ThisBuild / Test / parallelExecution := false

addCommandAlias(
  "coverageCoreFft",
  ";clean;coverage;coreJVM/test;fftJVM/test;lawsJVM/test;coverageAggregate"
)

// Gale: prefer an explicit override, then a sibling checkout, then a pinned git
// revision. Maven Central publication of gale-core is not assumed yet.
lazy val galeRevision = "d55fe2f97196a76ab7879e1a12f1e92403aeba06"
lazy val galeBuild: java.net.URI =
  sys.props
    .get("signal4s.gale.build")
    .map(path => file(path).getCanonicalFile.toURI)
    .getOrElse {
      val sibling = file("../gale").getCanonicalFile
      if (sibling.isDirectory) sibling.toURI
      else uri(s"https://github.com/canardlapin/gale.git#$galeRevision")
    }

lazy val galeCoreJVM = ProjectRef(galeBuild, "coreJVM")
lazy val galeCoreJS  = ProjectRef(galeBuild, "coreJS")

lazy val munitVersion = "1.3.0"
lazy val ujsonVersion = "4.4.3"

lazy val fixtureTestSupportSettings = Seq(
  Test / unmanagedSourceDirectories +=
    (ThisBuild / baseDirectory).value / "modules/test-support/jvm/src/main/scala",
  libraryDependencies += "com.lihaoyi" %% "ujson" % ujsonVersion % Test
)

lazy val jvmIncubatorVector = Seq(
  Compile / compile / javacOptions ++= Seq(
    "--add-modules",
    "jdk.incubator.vector",
    "--add-exports",
    "java.base/jdk.internal.vm.vector=ALL-UNNAMED"
  ),
  javaOptions ++= Seq("--add-modules", "jdk.incubator.vector"),
  Test / javaOptions ++= Seq(
    "--add-modules",
    "jdk.incubator.vector",
    "--enable-native-access=ALL-UNNAMED"
  ),
  run / javaOptions ++= Seq("--enable-native-access=ALL-UNNAMED"),
  Test / run / javaOptions ++= Seq("--enable-native-access=ALL-UNNAMED"),
  // Fork so module flags apply; keep repo root as cwd for fixture/receipt paths.
  Test / fork := true,
  run / fork := true,
  Test / run / fork := true,
  Test / baseDirectory := (ThisBuild / baseDirectory).value,
  Test / run / baseDirectory := (ThisBuild / baseDirectory).value
)

lazy val sharedSettings = Seq(
  Test / fork := false,
  libraryDependencies ++= Seq(
    "org.scalameta" %%% "munit"            % munitVersion % Test,
    "org.scalameta" %%% "munit-scalacheck" % munitVersion % Test
  )
)

lazy val jsSettings = Seq(
  scalaJSLinkerConfig ~= (_.withModuleKind(ModuleKind.CommonJSModule))
)

lazy val core =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/core"))
    .settings(sharedSettings)
    .settings(
      name := "signal4s-core",
      description :=
        "Sampled signals, kernels, and filter descriptions for Scala 3 on the JVM and Scala.js."
    )
    .jvmSettings(fixtureTestSupportSettings*)
    .jvmSettings(jvmIncubatorVector*)
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val coreJS  = core.js
lazy val coreJVM = core.jvm

lazy val fft =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/fft"))
    .dependsOn(core)
    .settings(sharedSettings)
    .settings(
      name := "signal4s-fft",
      description :=
        "Portable FFT plans, split-complex vectors, and real spectra for signal4s."
    )
    .jvmSettings(
      libraryDependencies += "com.github.wendykierp" % "JTransforms" % "3.1"
    )
    .jvmSettings(fixtureTestSupportSettings*)
    .jvmSettings(jvmIncubatorVector*)
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val fftJS  = fft.js
lazy val fftJVM = fft.jvm
lazy val design =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/design"))
    .dependsOn(core, fft)
    .settings(sharedSettings)
    .settings(
      name := "signal4s-design",
      description :=
        "FIR and IIR filter design (windowed-sinc, Butterworth, bilinear, ZPK/SOS)."
    )
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val designJS  = design.js
lazy val designJVM = design.jvm

lazy val laws =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/laws"))
    .dependsOn(core, fft, design)
    .settings(
      name := "signal4s-laws",
      description := "Reusable munit/ScalaCheck law bundles for signal4s.",
      Test / fork := false,
      libraryDependencies ++= Seq(
        "org.scalameta" %%% "munit"            % munitVersion,
        "org.scalameta" %%% "munit-scalacheck" % munitVersion
      )
    )
    .jsSettings(jsSettings)

lazy val lawsJS  = laws.js
lazy val lawsJVM = laws.jvm

// The coverage gate is scoped to the requested core + FFT surface. Laws still
// run against their instrumented dependencies, but their helper sources and
// design are not part of this metric.
designJVM / coverageEnabled := false
lawsJVM / coverageEnabled := false

lazy val ravel =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("modules/ravel"))
    .dependsOn(core)
    .settings(sharedSettings)
    .settings(
      name := "signal4s-ravel",
      description :=
        "Optional N-D / multi-channel adapters that apply 1-D signal4s ops along an axis."
    )
    .jvmConfigure(_.dependsOn(galeCoreJVM))
    .jsConfigure(_.dependsOn(galeCoreJS))
    .jsSettings(jsSettings)

lazy val ravelJS  = ravel.js
lazy val ravelJVM = ravel.jvm

// Public guides live in site-docs so the repository's audit and development
// notes under docs/ are not accidentally rendered or indexed.  The site is a
// local generation target only; publication is intentionally not configured.
lazy val docs =
  project
    .in(file("site"))
    .dependsOn(coreJVM, fftJVM, designJVM, ravelJVM)
    .settings(
      name := "signal4s-docs",
      publish / skip := true,
      mdocIn := (baseDirectory.value / ".." / "site-docs").getCanonicalFile,
      mdocExtraArguments := Seq(
        "--clean-target",
        "--check-link-hygiene",
        "--report-relative-paths"
      )
    )
    .enablePlugins(TypelevelSitePlugin)

lazy val root =
  project
    .in(file("."))
    .aggregate(
      coreJVM,
      coreJS,
      fftJVM,
      fftJS,
      designJVM,
      designJS,
      lawsJVM,
      lawsJS,
      ravelJVM,
      ravelJS
    )
    .settings(
      name := "signal4s",
      publish / skip := true,
      coverageMinimumStmtTotal := 85,
      coverageFailOnMinimum := true
    )
