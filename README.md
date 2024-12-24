# Black-Scholes Option Pricing Model

## Overview
This project is a Python-based application that calculates the price of European options using the Black-Scholes formula and Monte Carlo simulations. It integrates real-time market data and provides interactive visualizations to enhance understanding of options trading concepts.

## Features
- **Dynamic Pricing Models**: Calculate option prices using both the Black-Scholes analytical method and Monte Carlo simulations.
- **Real-Time Data Integration**: Fetch live stock data via YFinance and risk-free rates from the Federal Reserve Economic Data (FRED) API.
- **Interactive GUI**: Built with `tkinter`, the graphical interface allows users to visualize Profit/Loss (P/L) diagrams, Greeks, and adjust strike price, volatility, and time to maturity dynamically.
- **Detailed Visualizations**: Generate dynamic P/L diagrams with break-even analysis and probability of profit using `matplotlib`.

## Screenshots


## Key Features Explained

### Black-Scholes Analytical Pricing
The project calculates the price of American call and put options using the Black-Scholes formula, providing accurate results for standard market conditions.

### Monte Carlo Simulation
Monte Carlo simulations estimate the option price by simulating thousands of possible future stock price paths, offering a stochastic approach to valuation.

### Greeks Visualization
The application calculates and displays option Greeks, such as Delta, Gamma, Theta, Vega, and Rho, helping users understand the sensitivity of options to different parameters.

### Real-Time Data Integration
Stock prices and risk-free rates are fetched dynamically from YFinance and the FRED API, ensuring up-to-date calculations.

## Usage
1. Run the application
   ```
   python main.py
   ```
2.	Input the stock ticker, strike price, and time to maturity.
3.	Choose the option type (Call/Put) and view the P/L diagram, Greeks, and probability of profit.
