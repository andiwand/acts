add_unittest(CombinatorialKalmanFilter CombinatorialKalmanFilterTests.cpp)
add_unittest(TrackSelector TrackSelectorTests.cpp)
add_unittest(TrackParamsLookupAccumulator TrackParamsLookupAccumulatorTests.cpp)
add_unittest(RzMeasurementGrid RzMeasurementGridTests.cpp)
add_unittest(RzTransport RzTransportTests.cpp)
add_unittest(RzTrackFinder RzTrackFinderTests.cpp)
add_unittest(RzPropagation RzPropagationTests.cpp)
target_include_directories(
    ActsUnitTestRzPropagation
    PRIVATE ${PROJECT_SOURCE_DIR}/Core/src/TrackFinding/Rz
)
