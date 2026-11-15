function toggleNav() {
  document.querySelector('.nav-links').classList.toggle('nav-open');
}

function toggleTheme() {
  var current = document.documentElement.getAttribute('data-theme') || 'light';
  var next = current === 'dark' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', next);
  localStorage.setItem('theme', next);
}

document.addEventListener('DOMContentLoaded', function () {
  document.querySelector('.nav-toggle').addEventListener('click', toggleNav);
  document.getElementById('theme-toggle').addEventListener('click', toggleTheme);

  // Close nav when a link is clicked
  document.querySelectorAll('.nav-links a').forEach(function (link) {
    link.addEventListener('click', function () {
      document.querySelector('.nav-links').classList.remove('nav-open');
    });
  });
});
