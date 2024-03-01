---

title: "An Introduction to Monte-Carlo Simulations in Finance (with code implementation)"
date: 2024-02-27
categories: [The Monte-Carlo Series] 
tags: [simulations, derivative pricing, stochastic processes]
math: true

---

The Monte-carlo simulation is the workhorse of the financial engineering realm. It will not only hold it's ground against virtually every stochastic creature that you could throw at it, but sometimes might be your only available option. So... What really is a Monte-carlo Simulation?

> This post is only part 1 of a 2-part article introducing the Monte-Carlo simulation. You can check if the next part is available clicking [here](/categories/the-monte-carlo-series/).
{: .prompt-info }


## Monte-Carlo Methods Explained

The Monte-Carlo Methods are algorithms that rely on random sampling to aproximate numerical solutions to complex mathematical problems. The idea is that after many repeated random outcomes of an experiment, the set of results will tend on average to the real solution. With modern computers we can now easily simulate outrageous amounts of random paths of an event. Therefore, converting Monte-Carlo simulations into the most powerful tool in the arsenal of contemporary scientists (although machine-learning and AI may arrive to change that).

Generally the implementation of a Monte-Carlo goes like this: 

We define a model in which we want a dependent variable to be predicted, we specify the probability distribution of the indepent variables. And finally, we generate random values for the independent variables in each of the many simulations.


## Brief Overview of It's Origin

The modern version of the Monte-Carlo Method was invented by __Stanislaw Ulam__ while working in the Manhattan Project in Los Alamos National Laboratory –*Yes, the same guys who developed the atomic bomb* –. The story goes that while playing solitaire, Ulam realized it was easier to estimate the probability of a successful 52-card laid out by playing many games and writing down the results, than trying to calculate it with combinatorial mathematics. A few years down the road thanks to the growing capacity of calculation provided by the improvement of computers, this approach became feasible and the solution to many previously unsolvable problems.

# Monte-Carlo Simulations in Finance

By now you might be wondering: How does this seemingly unbeatable tool is used in finance?  
Well... let's find out.

Luckily for us, Monte-Carlo simulations can be used to solve problems within almost any financial model. Actually, It's very natural to reach for this tool, since in a dynamic setting most things in finance can be modeled as a stochastic process.

We'll now implement one of the simplest use cases, the pricing of a __european option__. We'll have to go over some math but I promise I'll keep it as simple as possible. And don't worry if you don't know what a stochastic process is. Stochastic is just a fancy word for random.


## Pricing a European Call with Monte-Carlo Integrations

### Setting up the problem

I will assume you have at least heard what an options is. In case you haven't, just think of it as the right to buy (or sell) an asset at a specific price and specific date. We call that specific price the __strike__, and that specific date the __maturity__. If you can buy the asset, that is a __call option__. If you can sell the asset, that is a __put option__.

The payoff of a european call option at maturity T reads:

$$
\begin{equation}
    C_T = max[S_T - K, 0]
\end{equation}
$$

That is, the maximum between 0 (the option is not exercised), or the difference between the strike and the spot price of the asset $$(S_t)$$ at maturity.

So let's say the company you're currently working for grants you some options to buy stock of said company a year from now, and of course, you want to know how much they are worth. How do you reasonably price this options? 

In order to price the option we can safely assume a condition of no arbitrage. This simply means that value of the option in the future should be equivalent to the price today. That is:

$$
\begin{equation}
    C_0 = \mathbb{E} \left[\frac{C_T}{e^{rT}}\right]
    \label{eq:nac}
\end{equation}
$$

* $$C_0$$ corresponds to value of the call today.
* $C_T$ is the value of the call at maturity. 
* __*T*__ is the maturity
* __*r*__ corresponds to the "risk-free" annual interest rate.
* $$e^{rT}$$ accounts for the __time value of money__.

> This last point is necessary because a sum of money is worth more in the present than in the future (due to the possibility of earning interest). The term $$\frac{C_T}{e^{rT}}$$ is called the discounted price of $$C_T$$.

If we knew what the value of $$C_T$$ is going to be, we could easily calculate the price today, just by discounting $$C_T$$. Nevertheless, as I don't have the ability to see the future –yet –and neither do you, we can only guess what the value $$C_T$$ will be. A good way to do this is by taking the expected value of the variable, hence the operator $$E[\cdot]$$.

In case you don't know what expected value means, just imagine it as the probability-weighted average of all possible outcomes of what's inside the brakets. The nature of this operator intuitively leads to think of Monte-Carlo simulations to try and solve it. So we just do that.

Do you remember the first thing we need to implement the algorithm?  
Exactly, a model that specifies the variable to be predicted and the underlying distributions.
The variable we want to estimate is simply the term of equation \eqref{eq:nac}. In the case of a N-path simulation, the estimator would be:
$$
\begin{equation}
    \bar{C}_0 = \frac{1}{N}\sum_{n=1}^{N}\left(\frac{C^{(n)}_{T}}{e^{rT}}\right)
    \label{eq:est}
\end{equation}
$$

For the underlying probability distribution we can make use of the famous Black-Scholes model, which specifies the dynamics under the risk-neutral measure of the market, or in this case of the underlying spot price $$(S_t)$$.

> If you're not familiar with the Black-Scholes model, don't freak out, we are not going to delve much into it, we'll just be using some of it's results. You don't really need to understand the model in full to follow along.

In the Black-Scholes model, the underlying returns follow a *Geometric Brownian Motion* or GBM. Which is nicely described by this stochastis differential equation:

$$
\begin{equation}
    dS_t = (r-q)S_t + \sigma S_t \cdot dW^Q_t
\end{equation}
$$

Fortunately for us, this SDE is solvable and yields a very convenient result for $$S_T$$:

$$
\begin{equation}
    S_T = S_0 e^{(r-q-\sigma^2/2)T +\sigma \sqrt{T}z}
\end{equation}
$$

+ *This result shows that $$\frac{S_T}{S_0}$$ is log-normally distributed. Here __z__ is a standard normal variable*

---

Now that we know our target to estimate and the distribution of the underlying, we can finally start with the fun part, the simulation itself. 

### Implementing the solution

In this section I will be writing code in Matlab, but feel free to use wathever your language of preference is.

> In case you want to get started with Matlab, Mathworks offers a 30-day free trial for all their tools in their [site](https://www.mathworks.com/).
{: .prompt-tip }

First of all, you might have notice that we need the values of certain parameters. Most of these can be observed in the market, in this instance I will provide you these values.

```matlab
clear
rng(34) %you can copy this seed to get identical results

So = 100;    %the spot price of the stock
r = 0.05;    %the risk-free rate
q = 0.03;    %the dividend the stock pays
sigma = 0.3; %the volatility of the underlying
```
*Notice that the volatility of the underlying is not observable. It is a model parameter and must be estimated as well. Let's assume we already did that work.*

We will also need to specify the details of the contract and the number of simulated paths:

```matlab
K = 120;    %the Strike price
T = 1;      %the Maturity in years
N = 10000;  %Number of paths
```
Now comes the simulation! For this, we will write a function to iterate over the number of paths and generate the corresponding random numbers. The function receives all the parameters we just defined and returns the discounted price estimate:

```matlab
function [call_estimate] = callMonteCarlo(So,r,q,K,T,sigma,N)

sum = 0;
```
Inside the function we defined an auxiliar variable called $$\texttt{sum}$$ that will accumulate the simulated discounted value of the call in each path, this variable is equivalent to the sum term in equation \eqref{eq:est}.

Next, we create a cycle where we generate a random sample of the underlying distribution as we use it to compute the corresponding payoff, then we just simply discount the result and we add it to $$\texttt{sum}$$.

```matlab
for i =1:N
    ST = So*exp((r-q -(sigma^2)/2 )*T + (sigma*sqrt(T)*randn())); 
    CT = max(ST-K,0); %we calculate the payoff
    YT = CT/exp(r*T); %we discount the payoff
    sum = sum + YT;
end
```
Last, we just divide the $$\texttt{sum}$$ by N –analog to eq. \eqref{eq:est} –and return the result.

```matlab
call_estimate = sum/N;
end
```
call this function in your __main.m__ and pass it all the previously defined parameters. displaying the result should yield:

```matlab
    5.8179
>>  
```
Eureka! now you know that today your options should be priced around __5.8__ each.

If you are not convinced yet, I got to tell you that there is actually a closed algebraic solution for this, given by the following formula:

$$
    C_0 = S_0 \cdot e^{-qT}\cdot \mathcal{N}(d_1) - K\cdot e^{-rT}\cdot \mathcal{N}(d_2)
$$

This is called the Black-Scholes formula for a call, and it's used to price european options. If we plug into it the same parameters we used in the simulation, the formula yields the price __5.8203__. Pretty close to our previous result!

* But wait, if there is a closed solution. Why did we bother doing simulations to price the option?

That's a valid point. The thing is, we went through a very simple example to get an idea of how to apply the Monte-Carlo method, and the fact there is a closed solution allow us to compare our results. Besides, this lovely formula only works with european vanilla options. Once we start working with slightly more complicated contracts, we don't have convenient formula to price them.

Hopefully you are now conviced this method works, and I'll be glad if this brief introduction was helpful to you. 

You can see the matlab code writen for this article in my [github](https://github.com/coffeeconomist).

> This post is only part 1 of a 2-part article introducing the Monte-Carlo simulation. You can check if the next part is available clicking [here](/categories/the-monte-carlo-series/).
{: .prompt-info }

