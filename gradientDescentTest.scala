import gradientDescent._

object gradientDescentTest {

  val precision = 1e-4

  val knownMinima = Map(
    "GD1" -> Point(Array(0.0, 0.0, 0.0)),
    "GD2" -> Point(Array(0.0, 0.0, 0.0)),
    "Boota" -> Point(Array(1.0, 3.0, 0.0)),
    "Matyas" -> Point(Array(0.0, 0.0, 0.0)),
    "HimmelblauMin1" -> Point(Array(3.0, 2.0, 0.0)),
    "HimmelblauMin2" -> Point(Array(-2.805118, 3.131312, 0.0)),
    "HimmelblauMin3" -> Point(Array(-3.779310, -3.283186, 0.0)),
    "HimmelblauMin4" -> Point(Array(3.584428, -1.848126, 0.0))
  )

  implicit def pToSet(p: Point): List[Point] = List(p)

  case class TestCase(
      name: String,
      function: Point => Double,
      startPoint: Point,
      expectedPoint: List[Point]
  )

  val testCases = List(
    TestCase(
      "Simple Quadratic: f(x, y) = (x² + y²)",
      testFunction.GD1,
      Point(Array(-9, 8)),
      knownMinima("GD1")
    ),
    TestCase(
      "Quadratic: f(x, y) = (3x²+y²-2xy)",
      testFunction.GD2,
      Point(Array(5, -5)),
      knownMinima("GD2")
    ),
    TestCase(
      "Booth Function: f(x, y) = (x + 2y -7)² + (2x + y - 5)²",
      testFunction.Boota,
      Point(Array(-9, 8)),
      knownMinima("Boota")
    ),
    TestCase(
      "Matyas Function: f(x, y) = 0.26(x² + y²) - 0.48xy",
      testFunction.Matyas,
      Point(Array(10, -10)),
      knownMinima("Matyas")
    ),
    TestCase(
      "Himmelblau's Function: f(x, y) = (x² + y - 11)² + (x + y² - 7)²",
      testFunction.Himmelblau,
      Point(Array(0, 0)),
      List(
        knownMinima("HimmelblauMin1"),
        knownMinima("HimmelblauMin2"),
        knownMinima("HimmelblauMin3"),
        knownMinima("HimmelblauMin4")
      )
    )
  )

  def formatNumber(value: Double): String = {
    if (s"%.5f".format(value).length > 8) s"%.5e".format(value)
    else s"%.5f".format(value)
  }

  def printTableHeader(): Unit = {
    println(
      "| %-20s | %-30s | %-15s | %-10s | %-12s |".format(
        "Method",
        "Found Point (x, y)",
        "Function Value",
        "Iterations",
        "Correctness"
      )
    )
    println(
      "|" + "-" * 22 + "|" + "-" * 32 + "|" + "-" * 17 + "|" + "-" * 12 + "|" + "-" * 14 + "|"
    )
  }

  def printTableRow(
      method: String,
      point: Point,
      value: Double,
      iterations: Int,
      correct: Boolean
  ): Unit = {
    val pointStr = s"(${formatNumber(point.x)}, ${formatNumber(point.y)})"
    val valueStr = formatNumber(value)
    val correctStr = if (correct) "YES" else "NO"

    println(
      "| %-20s | %-30s | %-15s | %-10d | %-12s |".format(
        method,
        pointStr,
        valueStr,
        iterations,
        correctStr
      )
    )
  }

  def runTests(): Unit = {
    println("=== Gradient Descent Method Testing ===")
    println(s"Comparison precision: $precision")
    println("=" * 50)

    testCases.foreach { testCase =>
      println(s"\nTesting function: ${testCase.name}")
      println(s"Start point: ${testCase.startPoint.printCords}")
      println("Expected minima:")
      testCase.expectedPoint.foreach(p => println(s"  ${p.printCords}"))
      println()

      val quadFunc = QuadFunc(testCase.function)

      val optimizationMethods = List(
        hConst,
        hSequence,
        hFunc,
        armijoRule(testCase.function),
        wolfeRule(testCase.function)
      )

      printTableHeader()

      optimizationMethods.foreach { scheduler =>
        val (minCords, iterations) =
          quadFunc.gradientDescent(testCase.startPoint, scheduler)

        val resultValue =
          Point(Array(minCords.x, minCords.y, testCase.function(minCords)))

        val isCorrect =
          testCase.expectedPoint.exists(p => p.equals(resultValue))

        printTableRow(
          scheduler.toString,
          minCords,
          resultValue.z,
          iterations,
          isCorrect
        )
      }
      println()
    }
  }

  def main(args: Array[String]): Unit = {
    runTests()

    println("\n=== Dichotomy Method Test ===")
    val dichotomyResult = dichotomy(fDichotomy, -1000, 1000)
    println(s"Found minimum of e^((x-10)²) at x = $dichotomyResult")
    println(s"Function value: ${fDichotomy(dichotomyResult)}")
    println(s"Expected: x ≈ 10.0, value = 1.0")
  }
}
