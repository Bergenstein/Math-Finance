> module Main where

> import Data.Complex
> import Data.List

> data OptionPricing = OptionPricing { rVal :: Double, s0Val :: Double, kVal :: Double, tVal :: Double }

> hestonCF :: OptionPricing -> Complex Double -> Double -> Double -> Double -> Double -> Double -> Complex Double
> hestonCF op u kappa_v theta_v sigma_v rho v0 =
>   let i = 0 :+ 1
>       c1 = kappa_v * theta_v
>       term = (rho * sigma_v * u * i - (kappa_v :+ 0))
>       c2 = negate (sqrt (term*term - sigma_v*sigma_v * ((-u*i) - u*u)))
>       c3 = ((kappa_v :+ 0) - rho * sigma_v * u * i + c2) / ((kappa_v :+ 0) - rho * sigma_v * u * i - c2)
>       h1 = (rVal op :+ 0) * u * i * (tVal op :+ 0) + (c1 / (sigma_v*sigma_v)) * (((kappa_v :+ 0) - rho * sigma_v * u * i + c2) * (tVal op :+ 0) - 2 * log ((1 - c3 * exp (c2 * (tVal op :+ 0))) / (1 - c3)))
>       h2 = (((kappa_v :+ 0) - rho * sigma_v * u * i + c2) / (sigma_v*sigma_v)) * ((1 - exp (c2 * (tVal op :+ 0))) / (1 - c3 * exp (c2 * (tVal op :+ 0))))
>   in exp (h1 + h2 * (v0 :+ 0))

> hestonIntegrationFunction :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double -> Double
> hestonIntegrationFunction op u kappa_v theta_v sigma_v rho v0 =
>   let i = 0 :+ 1
>       cf = hestonCF op (u - 0.5 * i) kappa_v theta_v sigma_v rho v0
>       integrand = exp (i * (u :+ 0) * log (s0Val op / kVal op)) * cf
>   in (1/(u*u+0.25)) * realPart integrand

> simpsonIntegrate :: (Double -> Double) -> Double -> Double -> Int -> Double
> simpsonIntegrate f a b n =
>   let h = (b - a) / fromIntegral n
>       sumTerms = f a + f b + sum [ (if even j then 2 else 4) * f (a + fromIntegral j * h) | j <- [1..n-1] ]
>   in (h / 3) * sumTerms

> hestonCallValue :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double
> hestonCallValue op kappa_v theta_v sigma_v rho v0 =
>   let f u = hestonIntegrationFunction op u kappa_v theta_v sigma_v rho v0
>       result = simpsonIntegrate f 0 100 2000
>   in max 0 (s0Val op - exp (-(rVal op)*tVal op) * sqrt (s0Val op * kVal op) / pi * result)

> mertonCFJump :: OptionPricing -> Complex Double -> Double -> Double -> Double -> Complex Double
> mertonCFJump op u lamb mu delta =
>   let i = 0 :+ 1
>       omega = -lamb * (exp (mu + 0.5*delta*delta) - 1)
>   in exp ((i*u*(omega :+ 0) + lamb*(exp (i*u*(mu :+ 0) - 0.5*u*u*delta*delta :+ 0) - 1)) * (tVal op :+ 0))

> mertonIntegrationFunction :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double
> mertonIntegrationFunction op u sigma lamb mu delta =
>   let i = 0 :+ 1
>       cf = mertonCFJump op (u - 0.5 * i) lamb mu delta
>       integrand = exp (i * (u :+ 0) * log (s0Val op / kVal op)) * cf
>   in (1/(u*u+0.25)) * realPart integrand

> mertonCallValue :: OptionPricing -> Double -> Double -> Double -> Double -> Double
> mertonCallValue op sigma lamb mu delta =
>   let f u = mertonIntegrationFunction op u sigma lamb mu delta
>       result = simpsonIntegrate f 0 100 2000
>   in max 0 (s0Val op - exp (-(rVal op)*tVal op) * sqrt (s0Val op * kVal op) / pi * result)

> batesCF :: OptionPricing -> Complex Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Complex Double
> batesCF op u kappa_v theta_v sigma_v rho v0 lamb mu delta =
>   hestonCF op u kappa_v theta_v sigma_v rho v0 * mertonCFJump op u lamb mu delta

> batesIntegrationFunction :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
> batesIntegrationFunction op u kappa_v theta_v sigma_v rho v0 lamb mu delta =
>   let i = 0 :+ 1
>       cf = batesCF op (u - 0.5 * i) kappa_v theta_v sigma_v rho v0 lamb mu delta
>       integrand = exp (i * (u :+ 0) * log (s0Val op / kVal op)) * cf
>   in (1/(u*u+0.25)) * realPart integrand

> batesCallValue :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
> batesCallValue op kappa_v theta_v sigma_v rho v0 lamb mu delta =
>   let f u = batesIntegrationFunction op u kappa_v theta_v sigma_v rho v0 lamb mu delta
>       result = simpsonIntegrate f 0 100 2000
>   in max 0 (s0Val op - exp (-(rVal op)*tVal op) * sqrt (s0Val op * kVal op) / pi * result)

> fft :: [Complex Double] -> [Complex Double]
> fft [x] = [x]
> fft xs =
>   let n = length xs
>       (evens, odds) = splitAt (n `div` 2) xs
>       evenFFT = fft evens
>       oddFFT = fft odds
>       twiddle k = exp (-2 * pi * (0 :+ 1) * fromIntegral k / fromIntegral n)
>       combine e o k = e + twiddle k * o
>   in zipWith combine evenFFT oddFFT [0..] ++ zipWith combine evenFFT oddFFT [n `div` 2 .. n-1]

> padToPowerOf2 :: [Complex Double] -> [Complex Double]
> padToPowerOf2 xs =
>   let n = length xs
>       nextPow2 = 2 ^ ceiling (logBase 2 (fromIntegral n))
>   in xs ++ replicate (nextPow2 - n) (0 :+ 0)

> batesCallFFT :: OptionPricing -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double -> Double
> batesCallFFT op kappa_v theta_v sigma_v rho v0 lamb mu delta =
>   let k = log (kVal op / s0Val op)
>       g = 1.0
>       n = floor (g * 4096)
>       eps = 1.0 / (g * 150.0)
>       eta = 2 * pi / (fromIntegral n * eps)
>       bVal = 0.5 * fromIntegral n * eps - k
>       us = [0..fromIntegral (n-1)]
>       vo = map (\j -> eta * j) us
>       i = 0 :+ 1
>       alpha = if s0Val op >= 0.95 * kVal op then 1.5 else 1.1
>       modcharFunc = if s0Val op >= 0.95 * kVal op
>                     then map (\v -> let vArr = v - (alpha+1) :+ 0 in exp (-(rVal op)*tVal op) * (batesCF op vArr kappa_v theta_v sigma_v rho v0 lamb mu delta / ((alpha*alpha+alpha) - v*v + i*(2*alpha+1)*v))) vo
>                     else let modcharFunc1 = map (\v -> let vArr = (v :+ 0) - (alpha :+ 0) - i in exp (-(rVal op)*tVal op) * (1/(1 + i*((v :+ 0) - (alpha :+ 0))) - exp ((rVal op)*tVal op)/(i*((v :+ 0) - (alpha :+ 0))) - batesCF op vArr kappa_v theta_v sigma_v rho v0 lamb mu delta/(((v :+ 0) - (alpha :+ 0))^2 - i*((v :+ 0) - (alpha :+ 0))))) vo
>                              modcharFunc2 = map (\v -> let vArr = (v :+ 0) + (alpha :+ 0) - i in exp (-(rVal op)*tVal op) * (1/(1 + i*((v :+ 0) + (alpha :+ 0))) - exp ((rVal op)*tVal op)/(i*((v :+ 0) + (alpha :+ 0))) - batesCF op vArr kappa_v theta_v sigma_v rho v0 lamb mu delta/(((v :+ 0) + (alpha :+ 0))^2 - i*((v :+ 0) + (alpha :+ 0))))) vo
>                          in zipWith (\m1 m2 -> (m1 - m2) * 0.5) modcharFunc1 modcharFunc2
>       simpsonW = 1 : [ (3 + (-1)^j) / 3 | j <- [2..n] ]
>       fftFuncList = zipWith3 (\v mc sw -> exp (i * (bVal :+ 0) * (v :+ 0)) * mc * (eta :+ 0) * (sw :+ 0)) vo modcharFunc simpsonW
>       fftFunc = padToPowerOf2 fftFuncList
>       payoff = fft fftFunc
>       callValM = if s0Val op >= 0.95 * kVal op
>                  then map (\z -> exp (-alpha*k)/pi * realPart z) payoff
>                  else map (\z -> realPart z / (sinh (alpha*k) * pi)) payoff
>       pos = floor ((k + bVal) / eps)
>   in (callValM !! pos) * s0Val op

> main :: IO ()
> main = do
>   let s0 = 100.0
>       k = 100.0
>       t = 1.0
>       r = 0.03
>       kappa_v = 1.5
>       theta_v = 0.02
>       sigma_v = 0.15
>       rho = 0.1
>       v0 = 0.01
>       lamb = 0.25
>       mu = -0.2
>       delta = 0.1
>       op = OptionPricing r s0 k t
>       batesPriceIntegration = batesCallValue op kappa_v theta_v sigma_v rho v0 lamb mu delta
>       batesPriceFFT = batesCallFFT op kappa_v theta_v sigma_v rho v0 lamb mu delta
>   putStrLn $ "Bates Call Price (Integration): " ++ show batesPriceIntegration
>   putStrLn $ "Bates Call Price (FFT): " ++ show batesPriceFFT