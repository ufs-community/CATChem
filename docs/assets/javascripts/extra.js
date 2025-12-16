// CATChem documentation JavaScript enhancements

// Add performance metrics animation
document.addEventListener('DOMContentLoaded', function() {
  // Animate performance metrics on scroll
  const performanceMetrics = document.querySelectorAll('.performance-metric');

  const observerOptions = {
    threshold: 0.3,
    rootMargin: '0px 0px -50px 0px'
  };

  const observer = new IntersectionObserver(function(entries) {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.style.opacity = '1';
        entry.target.style.transform = 'translateX(0)';
      }
    });
  }, observerOptions);

  performanceMetrics.forEach(metric => {
    metric.style.opacity = '0';
    metric.style.transform = 'translateX(-20px)';
    metric.style.transition = 'opacity 0.6s ease, transform 0.6s ease';
    observer.observe(metric);
  });
});

// Smooth scrolling for anchor links
document.addEventListener('click', function(e) {
  if (e.target.tagName === 'A' && e.target.getAttribute('href').startsWith('#')) {
    e.preventDefault();
    const targetId = e.target.getAttribute('href').substring(1);
    const targetElement = document.getElementById(targetId);
    if (targetElement) {
      targetElement.scrollIntoView({
        behavior: 'smooth',
        block: 'start'
      });
    }
  }
});

// CATChem Documentation Extra JavaScript
// No header branding - community project without institutional branding
