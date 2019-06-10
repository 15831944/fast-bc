#include <catch2/catch.hpp>

#include <brandes/VertexInfo.h>

using namespace fastbc::brandes;

TEST_CASE("Vertex info construtor", "[brandes]")
{
	VertexInfo<int, float> vi(4);

	REQUIRE_NOTHROW(vi.setBorderSPLength(3, 5.7f));

	REQUIRE_THROWS(vi.setBorderSPLength(4, 5.7f));

	SECTION("Copy constructor")
	{
		VertexInfo<int, double> viCopy(vi);

		REQUIRE(vi.getBorderSPLength(3) == viCopy.getBorderSPLength(3));

		viCopy.setBorderSPLength(3, 7.5);

		REQUIRE(vi.getBorderSPLength(3) != viCopy.getBorderSPLength(3));
	}

	SECTION("Copy assignment operator")
	{
		VertexInfo<int, double> viCopy = vi;

		REQUIRE(vi.getBorderSPLength(3) == viCopy.getBorderSPLength(3));

		viCopy.setBorderSPLength(3, 7.5);

		REQUIRE(vi.getBorderSPLength(3) != viCopy.getBorderSPLength(3));
	}
}

TEST_CASE("Getter and setters", "[brandes]")
{
	VertexInfo<short, float> vi(4);

	SECTION("Border SP length set/get")
	{
		REQUIRE(vi.getBorderSPLength(0) == 0.0f);

		vi.setBorderSPLength(0, 5.7f);

		REQUIRE(vi.getBorderSPLength(0) == 5.7f);
	}

	SECTION("Border SP count set/get")
	{
		REQUIRE(vi.getBorderSPCount(0) == 0);

		vi.setBorderSPCount(0, 7);

		REQUIRE(vi.getBorderSPCount(0) == 7);
	}

	SECTION("Minimum value getter")
	{
		vi.setBorderSPLength(0, 2.3f);
		vi.setBorderSPLength(1, 5.2f);
		vi.setBorderSPLength(2, 1.1f);
		vi.setBorderSPLength(3, 4.7f);

		REQUIRE(vi.getMinBorderSPLength() == 1.1f);
	}
}

TEST_CASE("Vertex info arithmetic operators", "[brandes]")
{
	VertexInfo<short, float> vi(4);
	vi.setBorderSPLength(0, 2.3f);
	vi.setBorderSPLength(1, 5.2f);
	vi.setBorderSPLength(2, 1.1f);
	vi.setBorderSPLength(3, 4.7f);
	vi.setBorderSPCount(0, 2);
	vi.setBorderSPCount(1, 5);
	vi.setBorderSPCount(2, 1);
	vi.setBorderSPCount(3, 4);

	SECTION("Normalize operator")
	{
		VertexInfo<short, float> viN(vi);
		viN.normalize();

		for (int i = 0; i < 4; i++)
		{
			REQUIRE(viN.getBorderSPLength(i) == (vi.getBorderSPLength(i) - vi.getMinBorderSPLength()));
			REQUIRE(viN.getBorderSPCount(i) == vi.getBorderSPCount(i));
		}
	}

	SECTION("Vextex info omogeneous operators")
	{
		VertexInfo<short, int> viD(4);
		viD.setBorderSPLength(0, 2);
		viD.setBorderSPLength(1, 5);
		viD.setBorderSPLength(2, 1);
		viD.setBorderSPLength(3, 4);
		viD.setBorderSPCount(0, 3);
		viD.setBorderSPCount(1, 4);
		viD.setBorderSPCount(2, 1);
		viD.setBorderSPCount(3, 4);

		SECTION("Squared distance operator")
		{
			REQUIRE(vi.squaredDistance(viD) == 2.62999964f);
		}

		SECTION("Sum")
		{
			VertexInfo<short, float> sum = vi + viD;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(sum.getBorderSPLength(i) == (vi.getBorderSPLength(i) + viD.getBorderSPLength(i)));
				REQUIRE(sum.getBorderSPCount(i) == (vi.getBorderSPCount(i) + viD.getBorderSPCount(i)));
			}

			vi += viD;

			REQUIRE(vi.getBorderSPLength(0) == 4.3f);
			REQUIRE(vi.getBorderSPLength(1) == 10.2f);
			REQUIRE(vi.getBorderSPLength(2) == 2.1f);
			REQUIRE(vi.getBorderSPLength(3) == 8.7f);
			REQUIRE(vi.getBorderSPCount(0) == 5);
			REQUIRE(vi.getBorderSPCount(1) == 9);
			REQUIRE(vi.getBorderSPCount(2) == 2);
			REQUIRE(vi.getBorderSPCount(3) == 8);
		}

		SECTION("Subtraction")
		{
			VertexInfo<short, float> sub = vi - viD;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(sub.getBorderSPLength(i) == (vi.getBorderSPLength(i) - viD.getBorderSPLength(i)));
				REQUIRE(sub.getBorderSPCount(i) == (vi.getBorderSPCount(i) - viD.getBorderSPCount(i)));
			}

			vi -= viD;

			REQUIRE(vi.getBorderSPLength(0) == (2.3f - (int)2));
			REQUIRE(vi.getBorderSPLength(1) == (5.2f - (int)5));
			REQUIRE(vi.getBorderSPLength(2) == (1.1f - (int)1));
			REQUIRE(vi.getBorderSPLength(3) == (4.7f - (int)4));
			REQUIRE(vi.getBorderSPCount(0) == -1);
			REQUIRE(vi.getBorderSPCount(1) == 1);
			REQUIRE(vi.getBorderSPCount(2) == 0);
			REQUIRE(vi.getBorderSPCount(3) == 0);
		}

		SECTION("Multiplication")
		{
			VertexInfo<short, float> mul = vi * viD;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(mul.getBorderSPLength(i) == (vi.getBorderSPLength(i) * viD.getBorderSPLength(i)));
				REQUIRE(mul.getBorderSPCount(i) == (vi.getBorderSPCount(i) * viD.getBorderSPCount(i)));
			}

			vi *= viD;

			REQUIRE(vi.getBorderSPLength(0) == 4.6f);
			REQUIRE(vi.getBorderSPLength(1) == 26.0f);
			REQUIRE(vi.getBorderSPLength(2) == 1.1f);
			REQUIRE(vi.getBorderSPLength(3) == 18.8f);
			REQUIRE(vi.getBorderSPCount(0) == 6);
			REQUIRE(vi.getBorderSPCount(1) == 20);
			REQUIRE(vi.getBorderSPCount(2) == 1);
			REQUIRE(vi.getBorderSPCount(3) == 16);
		}

		SECTION("Division")
		{
			VertexInfo<float, float> div(vi);
			div = div / viD;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(div.getBorderSPLength(i) == (vi.getBorderSPLength(i) / viD.getBorderSPLength(i)));
				REQUIRE(div.getBorderSPCount(i) == ((float)vi.getBorderSPCount(i) / viD.getBorderSPCount(i)));
			}

			vi /= viD;

			REQUIRE(vi.getBorderSPLength(0) == 1.15f);
			REQUIRE(vi.getBorderSPLength(1) == 1.04f);
			REQUIRE(vi.getBorderSPLength(2) == 1.1f);
			REQUIRE(vi.getBorderSPLength(3) == 1.175f);
			REQUIRE(vi.getBorderSPCount(0) == (short)(2 / 3));
			REQUIRE(vi.getBorderSPCount(1) == 1);
			REQUIRE(vi.getBorderSPCount(2) == 1);
			REQUIRE(vi.getBorderSPCount(3) == 1);
		}
	}

	SECTION("Vextex info eterogeneous operators")
	{
		float num = 3.4f;

		SECTION("Sum")
		{
			VertexInfo<short, float> sum = vi + num;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(sum.getBorderSPLength(i) == (vi.getBorderSPLength(i) + num));
				REQUIRE(sum.getBorderSPCount(i) == (short)(vi.getBorderSPCount(i) + num));
			}

			vi += num;

			REQUIRE(vi.getBorderSPLength(0) == (2.3f + num));
			REQUIRE(vi.getBorderSPLength(1) == (5.2f + num));
			REQUIRE(vi.getBorderSPLength(2) == (1.1f + num));
			REQUIRE(vi.getBorderSPLength(3) == (4.7f + num));
			REQUIRE(vi.getBorderSPCount(0) == (short)(2 + num));
			REQUIRE(vi.getBorderSPCount(1) == (short)(5 + num));
			REQUIRE(vi.getBorderSPCount(2) == (short)(1 + num));
			REQUIRE(vi.getBorderSPCount(3) == (short)(4 + num));
		}

		SECTION("Subtraction")
		{
			VertexInfo<short, float> sub = vi - num;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(sub.getBorderSPLength(i) == (vi.getBorderSPLength(i) - num));
				REQUIRE(sub.getBorderSPCount(i) == (short)(vi.getBorderSPCount(i) - num));
			}

			vi -= num;

			REQUIRE(vi.getBorderSPLength(0) == (2.3f - num));
			REQUIRE(vi.getBorderSPLength(1) == (5.2f - num));
			REQUIRE(vi.getBorderSPLength(2) == (1.1f - num));
			REQUIRE(vi.getBorderSPLength(3) == (4.7f - num));
			REQUIRE(vi.getBorderSPCount(0) == (short)(2 - num));
			REQUIRE(vi.getBorderSPCount(1) == (short)(5 - num));
			REQUIRE(vi.getBorderSPCount(2) == (short)(1 - num));
			REQUIRE(vi.getBorderSPCount(3) == (short)(4 - num));
		}

		SECTION("Multiplication")
		{
			VertexInfo<short, float> mul = vi * num;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(mul.getBorderSPLength(i) == (vi.getBorderSPLength(i) * num));
				REQUIRE(mul.getBorderSPCount(i) == (short)(vi.getBorderSPCount(i) * num));
			}

			vi *= num;

			REQUIRE(vi.getBorderSPLength(0) == (2.3f * num));
			REQUIRE(vi.getBorderSPLength(1) == (5.2f * num));
			REQUIRE(vi.getBorderSPLength(2) == (1.1f * num));
			REQUIRE(vi.getBorderSPLength(3) == (4.7f * num));
			REQUIRE(vi.getBorderSPCount(0) == (short)(2 * num));
			REQUIRE(vi.getBorderSPCount(1) == (short)(5 * num));
			REQUIRE(vi.getBorderSPCount(2) == (short)(1 * num));
			REQUIRE(vi.getBorderSPCount(3) == (short)(4 * num));
		}

		SECTION("Division")
		{
			VertexInfo<float, float> div(vi);
			div = div / num;

			for (int i = 0; i < 4; i++)
			{
				REQUIRE(div.getBorderSPLength(i) == (vi.getBorderSPLength(i) / num));
				REQUIRE(div.getBorderSPCount(i) == (vi.getBorderSPCount(i) / num));
			}

			vi /= num;

			REQUIRE(vi.getBorderSPLength(0) == (2.3f / num));
			REQUIRE(vi.getBorderSPLength(1) == (5.2f / num));
			REQUIRE(vi.getBorderSPLength(2) == (1.1f / num));
			REQUIRE(vi.getBorderSPLength(3) == (4.7f / num));
			REQUIRE(vi.getBorderSPCount(0) == (short)(2 / num));
			REQUIRE(vi.getBorderSPCount(1) == (short)(5 / num));
			REQUIRE(vi.getBorderSPCount(2) == (short)(1 / num));
			REQUIRE(vi.getBorderSPCount(3) == (short)(4 / num));
		}
	}
}

TEST_CASE("Vertex info compare operators", "[brandes]")
{
	VertexInfo<short, float> viA(4);
	viA.setBorderSPLength(0, 2.3f);
	viA.setBorderSPLength(1, 5.2f);
	viA.setBorderSPLength(2, 1.1f);
	viA.setBorderSPLength(3, 4.7f);
	viA.setBorderSPCount(0, 2);
	viA.setBorderSPCount(1, 5);
	viA.setBorderSPCount(2, 1);
	viA.setBorderSPCount(3, 4);

	VertexInfo<short, float> viAA(viA);

	VertexInfo<short, float> viB(4);
	viB.setBorderSPLength(0, 2.3f);
	viB.setBorderSPLength(1, 5.2f);
	viB.setBorderSPLength(2, 2.1f);
	viB.setBorderSPLength(3, 4.7f);
	viB.setBorderSPCount(0, 2);
	viB.setBorderSPCount(1, 5);
	viB.setBorderSPCount(2, 1);
	viB.setBorderSPCount(3, 4);

	VertexInfo<short, float> viC(4);
	viC.setBorderSPLength(0, 3.4f);
	viC.setBorderSPLength(1, 6.3f);
	viC.setBorderSPLength(2, 2.2f);
	viC.setBorderSPLength(3, 5.8f);
	viC.setBorderSPCount(0, 3);
	viC.setBorderSPCount(1, 6);
	viC.setBorderSPCount(2, 2);
	viC.setBorderSPCount(3, 5);

	SECTION("Equality comparator")
	{
		REQUIRE(viA == viA);
		REQUIRE(viA == viAA);
		REQUIRE_FALSE(viA == viB);
		REQUIRE_FALSE(viA == viC);

		REQUIRE(viA != viB);
		REQUIRE(viB != viC);
		REQUIRE_FALSE(viA != viAA);
	}

	SECTION("Inequality comparator")
	{
		REQUIRE(viA <= viA);
		REQUIRE(viA <= viB);
		REQUIRE(viA < viB);
		REQUIRE(viB < viC);
		REQUIRE(viA < viC);
		REQUIRE_FALSE(viB <= viA);
		REQUIRE_FALSE(viB < viA);

		REQUIRE(viA >= viAA);
		REQUIRE(viB >= viA);
		REQUIRE(viC > viB);
		REQUIRE(viB > viA);
		REQUIRE(viC > viA);
		REQUIRE_FALSE(viA > viAA);
		REQUIRE_FALSE(viA >= viB);
	}
}

TEST_CASE("Vertex info contribution distance", "[brandes]")
{
	VertexInfo<int, double> vi(3);
	vi.setBorderSPCount(0, 1);
	vi.setBorderSPLength(0, 10.0);
	vi.setBorderSPCount(1, 0);
	vi.setBorderSPLength(1, 0.0);
	vi.setBorderSPCount(2, 2);
	vi.setBorderSPLength(2, 6.0);

	VertexInfo<int, double> vib(3);
	vib.setBorderSPCount(0, 1);
	vib.setBorderSPLength(0, 10.0);
	vib.setBorderSPCount(1, 1);
	vib.setBorderSPLength(1, 5.0);
	vib.setBorderSPCount(2, 2);
	vib.setBorderSPLength(2, 6.0);

	REQUIRE(vi.contributionDistance(vib) == FASTBC_BRANDES_VERTEXINFO_PENALTY);
}