# first line: 95
    @classmethod
    def crawl_page(cls, url):
        crawler = LinkCrawler()
        file_urls = crawler.crawl(url)
        return file_urls
