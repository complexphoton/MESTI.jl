using Plots, LinearAlgebra, Statistics

function plot_and_compare_distribution(t)
    # Take the transmission matrix t, plot its transmission eigenvalue distribution, 
    # and compare it with the analytic distribution (bimodal DMPK distribution).

    # find transmission eigenvalue from svd of (adjoint(t)*t)
    _, tau, _ = svd(adjoint(t)*t)

    # define the edges of histogram
    bins = 0.0:0.02:1.0

    # initialize the counts array
    counts = zeros(length(bins) - 1)

    # count the number of in each bin
    for value in tau
        for i in 1:length(counts)
            if (value >= bins[i] && value < bins[i+1])
                counts[i] += 1
                break
            end
        end
    end
    counts[end] += sum(tau .== bins[end])

    # normalize to get the probability density (pdf)
    bin_width = step(bins)
    pdf = counts / (sum(counts) * bin_width)

    # get the midpoints of the bins for plotting
    bin_centers = (bins[1:end-1] .+ bins[2:end]) / 2

    # Plot the probability density distribution and compare it with the analytic distribution (bimodal DMPK distribution)
    p = bar(bin_centers, pdf, width=bin_width, label="Simulation", xlabel="Transmission eigenvalue", ylabel="Probability density",  size=(800, 600), tickfontsize=32, legendfontsize=32, guidefontsize=36, margins = 30Plots.px)
    plot!(p, bin_centers,mean(tau) ./ (2*bin_centers .* sqrt.(1 .- bin_centers)), label="Analytic result",  xlims = (0,1), linewidth=3, dpi= 300)
    display(p)
end