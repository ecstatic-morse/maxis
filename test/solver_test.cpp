#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>

#include "maxis/solver.hpp"

using namespace maxis;

class GeneticTest : public ::testing::Test {
public:
  GeneticTest() : wiki{0} {
    std::ifstream wiki_stream("test/wiki.mis");
    if (!wiki_stream) {
      throw std::runtime_error("Input graph not found");
    }
    wiki = parse_ascii_dimacs(wiki_stream);
  }

protected:
  Graph wiki;
};

TEST_F(GeneticTest, ValidIndependentSet) {
  EXPECT_TRUE(boost::edge(0, 1, wiki).second);
  EXPECT_TRUE(boost::edge(1, 0, wiki).second);
  for (int i = 0; i < 10; ++i) {
    auto m = Member::random(wiki);
    std::cout << m << std::endl;
    EXPECT_TRUE(m.is_valid(wiki));
  }
}

/*
TEST_F(GeneticTest, MaxisHeuristicMutator) {
    using namespace genetic;
    auto mut = MaxisHeuristicMutator{};
    mut.graph = &naaru;
    Phenotype pn{naaru.order()};
    std::fill(std::begin(*pn.chromosome), std::end(*pn.chromosome), 1);
    mut.mutate(pn);
    EXPECT_TRUE(naaru.is_independent_set(*pn.chromosome));

    mut.graph = &wiki;
    Phenotype pw{wiki.order()};
    std::fill(std::begin(*pw.chromosome), std::end(*pw.chromosome), 1);
    mut.mutate(pw);
    EXPECT_TRUE(wiki.is_independent_set(*pw.chromosome));
}
*/
