import numpy as np
from scipy.stats import norm
import yfinance as yf
import tkinter as tk
from tkinter import ttk
import requests
import os
from dotenv import load_dotenv

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
        "series_id": "DTB3",  # 3-Month Treasury Bill
        "api_key": FRED_API_KEY,
        "file_type": "json",
        "sort_order": "desc",  # Get the latest observation
        "limit": 1,
    }

    try:
        response = requests.get(FRED_URL, params=params)
        response.raise_for_status()  # Raise HTTPError for bad responses
        data = response.json()
        # Get the latest observation's value
        if "observations" in data and data["observations"]:
            rate = float(data["observations"][0]["value"]) / 100  # Convert from % to decimal
            return rate
    except Exception as e:
        print(f"Error fetching risk-free rate: {e}")
        return 0.05  # Default to 5% if API call fails

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
        ticker = ticker_entry.get()
        global current_stock_price
        current_stock_price, sigma = fetch_stock_info(ticker)

        # Update slider range dynamically based on stock price
        strike_price_slider.configure(from_=current_stock_price * 0.5,
                                       to=current_stock_price * 1.5,
                                       resolution=1)
        strike_price_slider.set(current_stock_price)  # Set default to current price

        result_text.set(f"Stock price retrieved: {current_stock_price:.2f}")
    except Exception as e:
        result_text.set(f"Error: {str(e)}")


def calculate_option():
    try:
        ticker = ticker_entry.get()
        K = float(strike_price_slider.get())  # Get strike price from slider
        T = float(maturity_entry.get()) / 365  # Convert days to years
        r = fetch_risk_free_rate()

        S, sigma = fetch_stock_info(ticker)

        model = BlackScholesModel(S, K, T, r, sigma)
        call_price = model.call_price()
        put_price = model.put_price()
        greeks = model.greeks()

        result_text.set(f"Ticker: {ticker}\n"
                        f"Current Stock Price: {S:.2f}\n"
                        f"Implied Volatility: {sigma:.2f}\n"
                        f"Call Option Price: {call_price:.2f}\n"
                        f"Put Option Price: {put_price:.2f}\n\n"
                        f"Greeks:\n" +
                        "\n".join(f"{key}: {value:.4f}" for key, value in greeks.items()))
    except Exception as e:
        result_text.set(f"Error: {str(e)}")


# Create the GUI
root = tk.Tk()
root.title("Black-Scholes Option Calculator")

# Input fields
ttk.Label(root, text="Stock Ticker:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
ticker_entry = ttk.Entry(root)
ticker_entry.grid(row=0, column=1, padx=5, pady=5)

# Button to fetch stock price and update slider
fetch_button = ttk.Button(root, text="Fetch Stock Info", command=update_strike_price_slider)
fetch_button.grid(row=0, column=2, padx=5, pady=5)

ttk.Label(root, text="Strike Price:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
strike_price_slider = tk.Scale(root, from_=50, to=150, orient=tk.HORIZONTAL, resolution=1)
strike_price_slider.grid(row=1, column=1, columnspan=2, sticky=tk.W + tk.E, padx=5, pady=5)

ttk.Label(root, text="Time to Maturity (days):").grid(row=2, column=0, sticky=tk.W, padx=5, pady=5)
maturity_entry = ttk.Entry(root)
maturity_entry.grid(row=2, column=1, padx=5, pady=5)

# Output field
result_text = tk.StringVar()
result_label = ttk.Label(root, textvariable=result_text, justify=tk.LEFT, anchor="w")
result_label.grid(row=6, column=0, columnspan=3, padx=5, pady=10)

# Calculate button
calculate_button = ttk.Button(root, text="Calculate", command=calculate_option)
calculate_button.grid(row=4, column=0, columnspan=3, pady=10)

# Run the GUI loop
root.mainloop()