for i in 1:5

    data = 10 * randn(440, 1) + 20 * sin(2 * pi * 0.01 * collect(1:440))
    cleaned = data .- Private.savitsky_golay(vec(data), 31, 2)

    @test mean(data[1:50]) > 6
    @test mean(cleaned[1:50]) < 1.5

    @test mean(data[100:150]) > 6
    @test mean(cleaned[100:150]) < 1.5

    # display(plot(data, lab = "raw"))
    # display(plot(data .- Private.savitsky_golay(vec(data), 11, 2), lab = "detrended", c = :red))
end

