> import System.Random
> import Control.Monad (replicateM)
> import Numeric.SpecFunctions (erf)
> import Data.List (intercalate)
> 
> data BlackScholes = BlackScholes { s0 :: Double, k :: Double, sigma :: Double, t :: Double, m :: Int, ite :: Int, r :: Double, dt :: Double }
> 
> newBlackScholes :: Double -> Double -> Double -> Double -> Int -> Int -> Double -> BlackScholes
> newBlackScholes s0' k' sigma' t' m' ite' r' = BlackScholes s0' k' sigma' t' m' ite' r' (t' / fromIntegral m')
> 
> normCDF :: Double -> Double
> normCDF x = 0.5 * (1 + erf (x / sqrt 2))
> 
> standardNormal :: IO Double
> standardNormal = do
>   u1 <- randomRIO (1e-10, 1)
>   u2 <- randomRIO (0, 1)
>   return (sqrt (-2 * log u1) * cos (2 * pi * u2))
> 
> simulatePaths :: BlackScholes -> IO [[Double]]
> simulatePaths bs = simulatePaths' (m bs) (replicate (ite bs) (s0 bs)) []
>   where
>     simulatePaths' 0 prev acc = return (reverse (prev:acc))
>     simulatePaths' n prev acc = do
>       zs <- replicateM (ite bs) standardNormal
>       let next = zipWith (\s z -> s * exp ((r bs - 0.5 * sigma bs^2) * dt bs + sigma bs * sqrt (dt bs) * z)) prev zs
>       simulatePaths' (n-1) next (prev:acc)
> 
> bsOptionMC :: BlackScholes -> Char -> IO Double
> bsOptionMC bs optype = do
>   paths <- simulatePaths bs
>   let lastPrices = last paths
>       payoff = case optype of
>         'C' -> map (\s -> max 0 (s - k bs)) lastPrices
>         'P' -> map (\s -> max 0 (k bs - s)) lastPrices
>         _   -> error "Invalid option type"
>       avg = sum payoff / fromIntegral (ite bs)
>   return (exp (-(r bs) * t bs) * avg)
> 
> bsOptionCF :: BlackScholes -> Char -> Double
> bsOptionCF bs optype =
>   let d1 = (log (s0 bs / k bs) + (r bs + 0.5 * sigma bs^2) * t bs) / (sigma bs * sqrt (t bs))
>       d2 = d1 - sigma bs * sqrt (t bs)
>   in case optype of
>     'C' -> s0 bs * normCDF d1 - k bs * exp (-(r bs) * t bs) * normCDF d2
>     'P' -> k bs * exp (-(r bs) * t bs) * normCDF (-d2) - s0 bs * normCDF (-d1)
>     _   -> error "Invalid option type"
> 
> printPaths :: BlackScholes -> Int -> IO ()
> printPaths bs n = do
>   if n <= 0 || n > ite bs then error "Invalid n" else return ()
>   paths <- simulatePaths bs
>   mapM_ (\(i, step) -> putStrLn ("Step " ++ show i ++ ": " ++ intercalate " " (map show (take n step)))) (zip [0..] paths)
> 
> main :: IO ()
> main = do
>   let s0' = 80.0
>       k' = 80.0
>       sigma' = 0.35
>       t' = 0.25
>       m' = 100
>       ite' = 10000
>       r' = 0.055
>       bs = newBlackScholes s0' k' sigma' t' m' ite' r'
>   callMC <- bsOptionMC bs 'C'
>   let callCF = bsOptionCF bs 'C'
>   putStrLn ("Call option price (Monte Carlo): " ++ show callMC)
>   putStrLn ("Call option price (Closed Form): " ++ show callCF)
