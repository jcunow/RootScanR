# img = load.image("C:/Users/jocu0013/Desktop/TamaraRoots/MinirhizotronR/binarized/mal_22_04_13_2b_1.png")
#
# img2 = imager::B(img)

#
# img_fft = fft(img2,inverse = T)
# # Compute magnitude spectrum
# magnitude <- Mod(img_fft)
# phase <- Arg(img_fft)
#
# phase_norm <- (phase + pi) / (2 * pi)
#
# # Display orientation map
# display(phase_norm, method = "raster", main = "Orientation Map")
#
# angle = atan2(Im(img_fft), Re(img_fft))
#
# # Calculate gradients using Sobel operator
# gradient <- imgradient(img,"xy",scheme = 3)
#
#
# y_grad = unlist(gradient[2],use.names = F)
# x_grad = unlist(gradient[1],use.names = F)
# # Extract gradient directions and magnitudes
# gradient_directions <- atan2(y_grad, x_grad)
# gradient_magnitudes <- sqrt(y_grad^2 + x_grad^2)
#
# # Visualize directionality lines
# # Example: Draw arrows based on gradient directions and magnitudes
# arrows_x <- cos(gradient_directions) * gradient_magnitudes
# arrows_y <- sin(gradient_directions) * gradient_magnitudes
#
# plot(img)
# arrows(x_grad, y_grad, x_grad + arrows_x, y_grad + arrows_y, length = 0.1)
