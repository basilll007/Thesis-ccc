/* ---------- Shared Dashboard UI Controller ---------- */
(function() {
  "use strict";

  // Theme management
  var root = document.documentElement;
  var themeBtn = document.getElementById('theme-toggle');
  var savedTheme = localStorage.getItem('ccc-theme');
  if (savedTheme) root.setAttribute('data-theme', savedTheme);

  function updateThemeIcon() {
    if (!themeBtn) return;
    var isDark = root.getAttribute('data-theme') === 'dark' ||
      (!root.getAttribute('data-theme') && window.matchMedia('(prefers-color-scheme: dark)').matches);
    themeBtn.textContent = isDark ? '☀' : '◐';
  }
  updateThemeIcon();

  if (themeBtn) {
    themeBtn.addEventListener('click', function() {
      var current = root.getAttribute('data-theme');
      var isDark = current === 'dark' || (!current && window.matchMedia('(prefers-color-scheme: dark)').matches);
      var next = isDark ? 'light' : 'dark';
      root.setAttribute('data-theme', next);
      localStorage.setItem('ccc-theme', next);
      updateThemeIcon();
      window.dispatchEvent(new CustomEvent('themechanged', { detail: { theme: next } }));
    });
  }

  // Scroll reveal
  var revealEls = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window && !window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
    var observer = new IntersectionObserver(function(entries) {
      entries.forEach(function(e) {
        if (e.isIntersecting) {
          e.target.classList.add('visible');
          observer.unobserve(e.target);
        }
      });
    }, { threshold: 0.05 });
    revealEls.forEach(function(el) { observer.observe(el); });
    setTimeout(function() { revealEls.forEach(function(el) { el.classList.add('visible'); }); }, 2000);
  } else {
    revealEls.forEach(function(el) { el.classList.add('visible'); });
  }

  // Global Lightbox
  window.openLightbox = function(src, caption) {
    var modal = document.getElementById('lightbox-modal');
    var img = document.getElementById('lightbox-img');
    var cap = document.getElementById('lightbox-caption');
    if (!modal || !img) return;
    img.src = src;
    if (cap) cap.textContent = caption || '';
    modal.classList.add('open');
    document.body.style.overflow = 'hidden';
  };

  window.closeLightbox = function() {
    var modal = document.getElementById('lightbox-modal');
    if (!modal) return;
    modal.classList.remove('open');
    document.body.style.overflow = '';
  };

  document.addEventListener('keydown', function(e) {
    if (e.key === 'Escape') window.closeLightbox();
  });
})();
