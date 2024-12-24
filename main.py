import os
import numpy as np
from scipy.stats import norm
import yfinance as yf
import tkinter as tk
from tkinter import ttk
import requests
from dotenv import load_dotenv
import matplotlib.pyplot as plt

load_dotenv()

class BlackScholesModel:
    def __init__(self, S, K, T, r, sigma):
        if S <= 0:
            raise ValueError("Stock price (S) must be greater than 0.")
        if K <= 0:
            raise ValueError("Strike price (K) must be greater than 0.")
        if T <= 0:
            raise ValueError("Time to maturity (T) must be greater than 0.")
        if sigma <= 0:
            raise ValueError("Volatility (sigma) must be greater than 0.")

        self.S = S
        self.K = K
        self.T = T
        self.r = r
        self.sigma = sigma
        self.d1 = self._calculate_d1()
        self.d2 = self.d1 - sigma * np.sqrt(T)

    def monte_carlo_simulation(self, num_simulations=10000, option_type="Call"):
        """
        Estimate option price using Monte Carlo simulation.
        :param num_simulations: Number of simulated price paths.
        :param option_type: "Call" or "Put".
        :return: Estimated option price.
        """

        np.random.seed(42)
        z = np.random.standard_normal(num_simulations)
        ST = self.S * np.exp((self.r - 0.5 * self.sigma**2) * self.T + self.sigma * np.sqrt(self.T) * z)

        # Payoff calculation
        if option_type == "Call":
            payoffs = np.maximum(ST - self.K, 0)
        elif option_type == "Put":
            payoffs = np.maximum(self.K - ST, 0)
        else:
            raise ValueError("Invalid option type. Use 'Call' or 'Put'.")

        option_price = np.exp(-self.r * self.T) * np.mean(payoffs)
        return option_price

    def _calculate_d1(self):
        return (np.log(self.S / self.K) + (self.r + 0.5 * self.sigma ** 2) * self.T) / (self.sigma * np.sqrt(self.T))

    def call_price(self):
        return self.S * norm.cdf(self.d1) - self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2)

    def put_price(self):
        return self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2) - self.S * norm.cdf(-self.d1)

    def greeks(self):
        delta_call = norm.cdf(self.d1)
        delta_put = delta_call - 1
        gamma = norm.pdf(self.d1) / (self.S * self.sigma * np.sqrt(self.T))
        vega = self.S * norm.pdf(self.d1) * np.sqrt(self.T) / 100
        theta_call = (-self.S * norm.pdf(self.d1) * self.sigma / (2 * np.sqrt(self.T))
                      - self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2)) / 365
        theta_put = (-self.S * norm.pdf(self.d1) * self.sigma / (2 * np.sqrt(self.T))
                     + self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2)) / 365
        rho_call = self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(self.d2) / 100
        rho_put = -self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(-self.d2) / 100

        return {
            "Delta (Call)": delta_call,
            "Delta (Put)": delta_put,
            "Gamma": gamma,
            "Vega": vega,
            "Theta (Call)": theta_call,
            "Theta (Put)": theta_put,
            "Rho (Call)": rho_call,
            "Rho (Put)": rho_put,
        }

def fetch_stock_info(ticker):
    stock = yf.Ticker(ticker)
    stock_data = stock.info
    S = stock_data.get('currentPrice', None)
    beta = stock_data.get('beta', 1.0)
    sigma = beta * 0.2

    if not S or S <= 0:
        raise ValueError("Invalid stock price retrieved.")
    return S, sigma

def update_strike_price_slider():
    try:
        ticker = ticker_entry.get().upper()
        current_stock_price, sigma = fetch_stock_info(ticker)

        strike_price_slider.config(from_=current_stock_price * 0.5,
                                   to=current_stock_price * 1.5,
                                   resolution=1)
        strike_price_slider.set(current_stock_price)

        result_text.set(f"Stock price retrieved: ${current_stock_price:.2f}")
    except Exception as e:
        result_text.set(f"Error: {str(e)}")

def fetch_risk_free_rate():
    """
    Fetch the 3-month Treasury yield from FRED.
    :return: Risk-free rate as a decimal (e.g., 0.05 for 5%)
    """
    FRED_API_KEY = os.getenv("FRED_API_KEY")
    if not FRED_API_KEY:
        raise ValueError("FRED_API_KEY not found. Please add it to your .env file.")
    FRED_URL = "https://api.stlouisfed.org/fred/series/observations"
    params = {
        "series_id": "DTB3",
        "api_key": FRED_API_KEY,
        "file_type": "json",
        "sort_order": "desc",
        "limit": 1,
    }

    try:
        response = requests.get(FRED_URL, params=params)
        response.raise_for_status()
        data = response.json()
        if "observations" in data and data["observations"]:
            rate = float(data["observations"][0]["value"]) / 100
            return rate
    except Exception as e:
        print(f"Error fetching risk-free rate: {e}")
        return 0.05



def calculate_and_plot_pl(option_type):
    try:
        ticker = ticker_entry.get().upper()
        K = float(strike_price_slider.get())
        T = float(maturity_entry.get()) / 365
        r = fetch_risk_free_rate()
        S, sigma = fetch_stock_info(ticker)

        model = BlackScholesModel(S, K, T, r, sigma)

        call_price = model.call_price()
        put_price = model.put_price()

        option_cost = call_price if option_type == "Call" else put_price
        breakeven_price = K + option_cost if option_type == "Call" else K - option_cost

        monte_carlo_call = model.monte_carlo_simulation(option_type="Call")
        monte_carlo_put = model.monte_carlo_simulation(option_type="Put")

        stock_prices = np.linspace(0.5 * S, 1.5 * S, 100)

        if option_type == "Call":
            pl = np.maximum(stock_prices - K, 0) - call_price
        else:
            pl = np.maximum(K - stock_prices, 0) - put_price

        plt.figure(figsize=(12, 7))
        plt.plot(stock_prices, pl, label=f'{option_type} P/L', color='blue', linewidth=2)

        plt.fill_between(stock_prices, pl, where=(pl > 0), interpolate=True, color='green', alpha=0.3, label='Profit Region')
        plt.fill_between(stock_prices, pl, where=(pl <= 0), interpolate=True, color='red', alpha=0.3, label='Loss Region')

        plt.axhline(0, color='black', linestyle='--', label='Break-Even Line')
        plt.axvline(breakeven_price, color='orange', linestyle='--', label=f'Break-Even Price (${breakeven_price:.2f})')

        plt.title(f'{option_type} Option P/L Diagram', fontsize=16)
        plt.xlabel(f'{ticker.upper()} Stock Price at Expiration', fontsize=14)
        plt.ylabel('Profit / Loss ($)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend(fontsize=12)

        POP = 1 - norm.cdf(model.d2) if option_type == "Call" else norm.cdf(-model.d2)
        plt.text(0.7 * stock_prices.max(), 0.7 * pl.max(), f"Probability of Profit: {POP:.2%}", fontsize=12, color='black',
                 bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))

        plt.show()

    except Exception as e:
        result_text.set(f"Error: {str(e)}")


def update_plot_and_greeks(*args):
    """
    Recalculate and update the P/L plot and Greeks dynamically based on sliders.
    """
    try:
        ticker = ticker_entry.get().upper()
        current_stock_price, _ = fetch_stock_info(ticker)
        K = float(strike_price_slider.get())
        T = float(maturity_slider.get()) / 365
        r = fetch_risk_free_rate()
        S = float(stock_price_slider.get())
        sigma = float(volatility_slider.get()) / 100

        model = BlackScholesModel(S, K, T, r, sigma)
        call_price = model.call_price()
        put_price = model.put_price()
        greeks = model.greeks()

        option_cost = call_price if option_type.get() == "Call" else put_price
        breakeven_price = K + option_cost if option_type.get() == "Call" else K - option_cost

        POP = 1 - norm.cdf(model.d2) if option_type.get() == "Call" else norm.cdf(-model.d2)

        stock_prices = np.linspace(0.5 * current_stock_price, 1.5 * current_stock_price, 100)

        if option_type.get() == "Call":
            pl = np.maximum(stock_prices - K, 0) - call_price
            label = f'Call P/L (K={K})'
        else:
            pl = np.maximum(K - stock_prices, 0) - put_price
            label = f'Put P/L (K={K})'

        ax.clear()
        ax.plot(stock_prices, pl, label=label, color='blue')
        ax.axhline(0, color='black', linestyle='--', linewidth=1, label='Break-Even Line')
        ax.axvline(breakeven_price, color='green', linestyle='--', linewidth=1, label=f'Break-Even (${breakeven_price:.2f})')
        ax.set_xlabel(f'{ticker.upper()} Stock Price at Expiration')
        ax.set_ylabel('Profit / Loss')
        ax.set_title(f'{option_type.get()} Option P/L Diagram')
        ax.legend()
        ax.grid()
        canvas.draw()

        greeks_display = "\n".join([f"{key}: {value:.4f}" for key, value in greeks.items()])
        greeks_text.set(f"{option_type.get()} Greeks:\n{greeks_display}\n"
                        f"Cost of Option: ${option_cost:.2f}\n"
                        f"Break-Even Price: ${breakeven_price:.2f}\n"
                        f"Probability of Profit: {POP:.2%}")

    except Exception as e:
        result_text.set(f"Error: {str(e)}")

root = tk.Tk()
root.title("Black-Scholes Option Calculator")

ttk.Label(root, text="Stock Ticker:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
ticker_entry = ttk.Entry(root)
ticker_entry.grid(row=0, column=1, padx=5, pady=5)

fetch_button = ttk.Button(root, text="Fetch Stock Info", command=update_strike_price_slider)
fetch_button.grid(row=0, column=2, padx=5, pady=5)

ttk.Label(root, text="Strike Price:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
strike_price_slider = tk.Scale(root, from_=50, to=150, orient=tk.HORIZONTAL, resolution=1)
strike_price_slider.grid(row=1, column=1, columnspan=2, sticky=tk.W + tk.E, padx=5, pady=5)

ttk.Label(root, text="Time to Maturity (days):").grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
maturity_entry = ttk.Entry(root)
maturity_entry.grid(row=2, column=1, padx=5, pady=5)

ttk.Label(root, text="Option Type:").grid(row=3, column=0, sticky=tk.W, padx=5, pady=5)
option_type = tk.StringVar(value="Call")
option_type_dropdown = ttk.OptionMenu(root, option_type, "Call", "Call", "Put")
option_type_dropdown.grid(row=3, column=1, padx=5, pady=5)

ttk.Label(root, text="Greeks & Prices:").grid(row=4, column=0, sticky=tk.W, padx=5, pady=5)
greeks_text = tk.StringVar()
greeks_label = ttk.Label(root, textvariable=greeks_text, justify=tk.LEFT, anchor="w")
greeks_label.grid(row=4, column=1, columnspan=2, padx=5, pady=5)

plot_button = ttk.Button(root, text="Plot P/L", command=lambda: calculate_and_plot_pl(option_type.get()))
plot_button.grid(row=7, column=0, columnspan=3, pady=10)

result_text = tk.StringVar()
result_label = ttk.Label(root, textvariable=result_text, justify=tk.LEFT, anchor="w")
result_label.grid(row=6, column=0, columnspan=3, padx=5, pady=10)

root.mainloop()