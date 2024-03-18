---

title: "An Introduction to Monte-Carlo Simulations in Finance #2: Measure of Accuracy"
date: 2024-03-17
categories: [The Monte-Carlo Series] 
tags: [simulations, derivative pricing, stochastic processes, tutorials, matlab]
math: true

---

Previously we introduced one of the simplest possible implementations of a Monte-Carlo simulation to price a european option. It is –as we've seen –certainly one of the most powerful tools you can learn in finance. However, it doesn't come without it's shortcomings. In this brief post we will be taking a look at one of the aspects that makes Monte-Carlo simulations not so great all the time. They are slow.

If you haven't already, you can read first part of this article [here](/posts/montecarlo-simulations-1). I highly encourage you to check it out since we will be modyfing the code written for part 1.

## Error of the estimation

You might be wondering what do I mean by slow. Well... in order to understand the method's speed issue, we firstly need to introduce the concept of error and accuracy of the estimation.

Let's take a look at the results from our previous implementation:

```matlab
    5.8179 % Monte-Carlo estimated price
>>  
    5.8203 % Black-Scholes formula price
>>
```
The estimated price and the price solution are quite close in this example. We can define the error of the estimation as the absolute value of the difference between these two values:

$$
\begin{equation*}
    \epsilon =  \left| \bar{C}_0^{MC} - C_0^{BS} \right| = | 5.8179 - 5.8203| = 0.0024 \\ 
\end{equation*}
$$

So the error of this estimation would be __0.0024__, very small indeed. Nevertheless, see what happens when we try estimating the price with different random seeds:

```matlab
rng(34)
>>  0.0024

rng(150)
>>  0.1341

rng(27)
>>  0.0259
```
The error wildly varies depending on the seed we use. This is concerning since we don't always have access to a closed formula to compare and calculate how big the error is, or you might not be able to choose a seed (for instance, when implementing parallelization). Plus, in finance the difference in price of a few basis points might mean a lot of money. In summary, the error is not a good way of assessing the reliability of our estimation.

Then... How can we trust our price estimator? For this reason, we need a measure of precision of the estimation.

## Understanding Accuracy

A good metric for accuracy is the estimator's standard deviation. Remember that our Monte-Carlo estimator reads:

$$
\begin{equation}
    \bar{C}_0 = \frac{1}{N}\sum_{n=1}^{N}\left(\frac{C^{(n)}_{T}}{e^{rT}}\right)
    \label{eq:est}
\end{equation}
$$

First, let's check if our estimator is biased. For this, we calculate the expectancy:

$$
\begin{equation*}
    \mathbb{E}\left[ \bar{C}_0 \right] = \frac{1}{N} \sum_{n=1}^{N} \mathbb{E}\left( \frac{C^{(n)}_{T}}{e^{rT}} \right)
\end{equation*}
$$

Since the discounted payoffs are identically distributed:

$$
\begin{equation}
    \mathbb{E}\left[ \bar{C}_0 \right] = \frac{1}{N} \sum_{n=1}^{N} \mathbb{E}\left( \frac{C_{T}}{e^{rT}} \right) = \mathbb{E}^Q\left( \frac{C_{T}}{e^{rT}} \right)= C_0
\end{equation}
$$

This means that our estimator is unbiased!  
Now, let's check the variance:

$$
\begin{equation*}
    Var\left[ \bar{C}_0 \right] = Var\left[ \frac{1}{N}\sum_{n=1}^{N}\left(\frac{C^{(n)}_{T}}{e^{rT}}\right)\right] = \frac{1}{N^2} \sum_{n=1}^N Var\left[\frac{C^{(n)}_{T}}{e^{rT}} \right]
\end{equation*}
$$

>Notice that unlike the first moment, the variance does depend on the number of paths

Again, since the payoffs are indentically distributed:

$$
\begin{equation}
    Var\left[ \bar{C}_0 \right] = \frac{1}{N^2} \sum_{n=1}^N Var\left[\frac{C_{T}}{e^{rT}} \right] = \frac{C_{T}}{N\cdot e^{rT}} 
\end{equation}
$$

Therefore, the measure of accuracy reads:

$$
\begin{equation}
    \Delta \bar{C}_0 = \sqrt{Var[\bar{C}_0]} = \frac{1}{\sqrt{N}} Std\left[\frac{C_{T}}{e^{rT}} \right]
    \label{eq:acc}
\end{equation}
$$

We can draw a few conclusions from this result. Notice that accuracy increases as $$\Delta \bar{C}_0$$ decreases. More so, the convergence rate is $$ \frac{1}{\sqrt{N}}$$.

## The speed issue

In order to have a more reliable price we can target a certain precision $$\Delta\bar{C}_0$$. But decreasing this metric comes at a cost. As we've seen, the convergence rate is $$ \frac{1}{\sqrt{N}}$$, this means that we need to exponentially increase the number of simulated paths if we want a linear increase in precision. In practice, this is a problem due to the fact that computational resources are finite, we can only have so much memory and procesing power. In other words, the problem is not the number of paths, but the computational time that increases linearly with the number of paths. Hence, a precise Monte-Carlo simulation will be slow.

The following figure illustrates the evolution of error within a band of plus-minus one measure of accuracy as the number of paths increments.

![accuracy-convergence](/assets/img/posts/montecarlo2/accuracy-convergence.png)

Hopefully, you can see that when approaching a smaller error, we need to exponentially increase the number of simulated paths in order to obtain a small gain in precision. And as we now, the more paths, the slower the simulation.

### Implementation

Let's see this effect in practice with code. We can modify the function from our previous implementation that we used to compute the price of a call:

```matlab
function [call_estimate, delta] = callMonteCarlo(So,r,q,K,T,sigma,N)
tic
sum = 0; %accumulador of payoffs
sum2 = 0; %accumulator of squared payoffs

for i =1:N
    ST = So*exp((r-q -(sigma^2)/2 )*T + (sigma*sqrt(T)*randn())); 
    CT = max(ST-K,0); %we calculate the payoff
    YT = CT/exp(r*T); %we discount the payoff
    sum = sum + YT;
    sum2 = sum2 + YT^2;
end
call_estimate = sum/N;     %Price estimator
var = sum2/N - (sum/N)^2;  %variance estimator
delta = sqrt(var)/sqrt(N); %measure of accuracy
toc
end
```

> You can see the entire code [here](https://github.com/coffeeconomist/tutorials-code/tree/main/montecarlo-1-2).
{: .prompt-info }

As you can see, we added a new accumulator for the squared discounted payoff, we then use it to estimate the variance of said random variable. Finally, we can simply use the expression from equation \eqref{eq:acc} to obtain our measure of accuracy (delta).

> In Matlab you can use "tic" and "toc" to display the elapsed time of execution for a block of code.
{: .prompt-tip }

If we display the error, accuracy and elapsed time for the same seeds we used before, we obtain the following results:

| Seed        | Error       | $$\Delta\bar{C}_0$$ | Ex. Time (milliseconds) |
| ----------- | ----------- | ----------- |         -----------  |
| rng(34)     | 0.0024      | 0.1519      |        1.429 |  
| rng(150)    | 0.1341      | 0.1474      |        0.845 |
| rng(27)     | 0.0259      | 0.1510      |        0.885 |

As you can see, the measure of accuracy $$\Delta\bar{C}_0$$ results in a very stable indicator of the price estimation reliability.


Now, let's see how the output changes when we increase the number of paths from 10,000 to 1 million:

| Seed        | Error       | $$\Delta\bar{C}_0$$ | Ex. Time (milliseconds) |
| ----------- | ----------- | ----------- |         -----------  |
| rng(34)     | 5.6077e-04      | 0.0150      |    89.913 | 
| rng(150)    | 0.0220     | 0.0150      |    109.253 |
| rng(27)     | 0.0025     | 0.0151      |    86.625|


As expected, for the precision to drop by __one__ significant digit we had to multiply the number of paths by __100__. However, doing this dramatically increased the execution time for our estimation. It's easy to imagine how this trade-off could be an issue with more complex and demanding implementations.

## Conclusions

You could think that 100 milliseconds aren't a big deal and you would be right. But the point is, once you start introducing more complexity to the implemented models, the calculation time may increase manyfold. For instance, you might have to use a local volatility model, or have many correlated underlying assets. This would not only impact on the calculations per path, but also on the memory usage. Additionally, targeting a very small error band might end up turning those initial 100 miliseconds into several minutes, or even hours. 

For all these reasons, when implementing Monte-Carlo simulations it's important to have a focus on writing efficient code and implementing speed-up techniques.

> There is one simple trick you can use to speed up a simulation like this in Matlab, by using it's parallel computing toolkit. You can run cycles in parallel by replacing all __for__ statements by __parfor__. As an activity, try using both versions and compare the elapsed time for the simulations.
{: .prompt-tip }

If case you found this post interesting, I will be writing more on this topic in the future. You can get notify whenever I publish a new entry of this blog by following me on [linkedin](https://www.linkedin.com/in/gdelpinocastro/).