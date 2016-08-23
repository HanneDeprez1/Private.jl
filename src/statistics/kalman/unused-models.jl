########################
#
# These are models that I tried during the development of the kalman filter
# In the end I settled on the basic model included in kalman.jl
#
########################


# function EASSR_model(f::Number, fs::Number; x0::Vector=[0.0, 0, 0, 150, 0, 0, 0, 0, 0, 0, 1], W::Number=1.0)
#
#     Δt = 1 / fs
#     ω = 2 * pi * f
#
#     ## Process transition and noise covariance
#
#     a = zeros(11,11)
#     a[1:2, 1:2] = eye(2)
#     a[4:11, 3:10] = eye(8)
#     a[3, 11] = 1
#
#     F = a
#
#     V = diagm([1e-10, 1e-10, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001])  # We expect little variation in the state matrix
#
#     ## Observation and noise covariance
#
#     function G1(k);  cos(ω * Δt * k); end
#     function G2(k); -sin(ω * Δt * k); end
#     G = reshape([G1, G2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1, 11)
#     W = diagm([W])   # Noise Covariance
#
#     ## Inital guesses at state and error covariance
#
#     P0 = diagm([0.2, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
#
#     ASSR_Model = StateSpaceModel(F, V, G, W, x0, P0)
#
# end

# function click_model(f::Number, fa::Number, fs::Number; x0::Vector=vcat(vec(0.1 * ones(6, 1)), vec(zeros(4, 1)), vec(20* ones(10, 1))),
#     W::Number=0.0001, V = diagm(vcat(vec(1e-10* ones(10, 1)), vec(0.001* ones(10, 1)))))
#
#     Δt = 1 / fs
#     ω = 2 * pi * f
#     ωa = 2 * pi * fa
#
#     ## Process transition and noise covariance
#
#     F = eye(20)
#     # V = diagm(vcat(vec(1e-10* ones(10, 1)), vec(0.001* ones(10, 1))))
#
#     ## Observation and noise covariance
#
#     names = (:G1, :G3, :G5, :G7, :G9)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             cos($freq * $ωa * $Δt * k)
#         end
#     end
#
#     names = (:G2, :G4, :G6, :G8, :G10)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             -sin($freq * $ωa * $Δt * k)
#         end
#     end
#
#     names = (:M1, :M3, :M5, :M7, :M9)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             cos($freq * $ω * $Δt * k)
#         end
#     end
#
#     names = (:M2, :M4, :M6, :M8, :M10)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             -sin($freq * $ω * $Δt * k)
#         end
#     end
#
#     G = reshape([M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10], 1, 20)
#     W = diagm([W])   # Noise Covariance
#
#     ## Inital guesses at state and error covariance
#
#     P0 = 20 * eye(20)
#
#     ASSR_Model = StateSpaceModel(F, V, G, W, x0, P0)
#
# end



# function ftest_model(f::Number, num_noise::Int, noise_freq::Number, fs::Number, W::Number;
#     x0::Vector = zeros(2 + (2 * num_noise), 1),  V = diagm(1e-10 * ones(2 + (2 * num_noise), 1)))
#
#     Δt = 1 / fs
#     ω = 2 * pi * f
#     ωa = 2 * pi * fa
#
#     ## Process transition and noise covariance
#
#     F = eye(20)
#     # V = diagm(vcat(vec(1e-10* ones(10, 1)), vec(0.001* ones(10, 1))))
#
#     ## Observation and noise covariance
#
#     names = (:G1, :G3, :G5, :G7, :G9)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             cos($freq * $ωa * $Δt * k)
#         end
#     end
#
#     names = (:G2, :G4, :G6, :G8, :G10)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             -sin($freq * $ωa * $Δt * k)
#         end
#     end
#
#     names = (:M1, :M3, :M5, :M7, :M9)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             cos($freq * $ω * $Δt * k)
#         end
#     end
#
#     names = (:M2, :M4, :M6, :M8, :M10)
#     freqs = (1,   2,   3,   4,   5)
#     for (name, freq) in zip(names, freqs)
#         @eval function ($name)(k)
#             -sin($freq * $ω * $Δt * k)
#         end
#     end
#
#     G = reshape([M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10], 1, 20)
#     W = diagm([W])   # Noise Covariance
#
#     ## Inital guesses at state and error covariance
#
#     P0 = 20 * eye(20)
#
#     ASSR_Model = StateSpaceModel(F, V, G, W, x0, P0)
#
# end
