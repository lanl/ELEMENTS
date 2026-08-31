#include <gtest/gtest.h>
#include <matar.h>

int main(int argc, char** argv)
{
    MATAR_INITIALIZE(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MATAR_FINALIZE();
    return result;
}
