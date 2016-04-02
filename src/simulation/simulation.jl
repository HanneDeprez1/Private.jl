using EEG, Private, JLD, ProgressMeter, ParallelAccelerator

ldf = load("/Users/rluke/Documents/Transfer-BTSync/Headmodels/Leadfield-Bilateral-CR20.jld")["ldf"]
a = read_SSR("/Users/rluke/Data/EEG/NH/Example/Example-40Hz.bdf")

l = size(a.data, 1)

ldf.L = permutedims(ldf.L, [3, 2 ,1])

result = zeros(l, 64)

println(size(result))

@acc @showprogress for i in 1:size(ldf.L, 1)
    for j in 1:3

        sig = 20 * ((rand(l, 1) * 2) - 1)

        proj = ldf.L[:, j, i]

        result += sig * proj'

    end
end

result /= size(ldf.L, 1) * 3

a.data = result;

a = rereference(a, "Cz")
a = extract_epochs(a)
a = epoch_rejection(a)
a = hotelling(a)
