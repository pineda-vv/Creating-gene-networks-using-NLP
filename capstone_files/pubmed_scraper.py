from selenium import webdriver
import selenium
import time
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
import numpy as np
import string
import pandas as pd
from collections import defaultdict
from pymongo import MongoClient


def exploratory_search(link):
    browser = webdriver.Firefox()
    browser.get(link)
    search_box = browser.find_element_by_css_selector(
        "#maincontent > div > div:nth-child(5) > div:nth-child(1) > div.rslt > p > a"
    )
    search_box.click()
    # search_box.send_keys(zipcode) #works up to here
    # search_button = browser.find_element_by_css_selector(
    #         "body > div.section.section-green > div > div.side-image-content > form > input.btn-white"
    #     )
    # search_button.click()
    # time.sleep(3)
    # search_button2 = browser.find_element_by_css_selector(
    #         "body > div:nth-child(3) > div > ul > li:nth-child(1) > a"
    #     )
    # search_button2.click()
    # time.sleep(3)   # Works to here now

    # salary = browser.find_element_by_xpath(
    #         "/html/body/div[2]/div[1]/div/table/tbody/tr[2]"
    #     )
    #
    # income = str(salary.text).split()
    # browser.close()
    # return income


if __name__ == "__main__":
    exploratory_search("https://www.ncbi.nlm.nih.gov/pubmed/?term=%22PLOS+Genetics%22")
